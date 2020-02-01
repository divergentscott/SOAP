#include "bdfIO.h"
#include <vtkCellType.h>
#include <algorithm>
#include <fstream>
#include <iostream>
#include <sstream>

namespace d3d {

namespace io {
namespace {
const int defaultPID = 1;

const auto triName = "CTRIA3  ";
const auto quadName = "CQUAD4  ";
const auto tetraName = "CTETRA  ";
const auto hexaName = "CHEXA   ";
const auto pyraName = "CPYRA   ";

const CellElement hexaElem =
    CellElement("CHEXA", VTK_HEXAHEDRON, 8, 3, CellType::CHEXA);
const CellElement tetraElem =
    CellElement("CTETRA", VTK_TETRA, 4, 3, CellType::CTETRA);
const CellElement quadElem =
    CellElement("CQUAD4", VTK_QUAD, 4, 2, CellType::CQUAD4);
const CellElement triElem =
    CellElement("CTRIA3", VTK_TRIANGLE, 3, 2, CellType::CTRIA3);
const CellElement pyraElem =
    CellElement("CPYRA", VTK_PYRAMID, 5, 3, CellType::CPYRA);

const std::vector<CellElement> allCellElements{hexaElem, tetraElem, quadElem,
                                               triElem, pyraElem};

const int standardCharSpace = 8;
const int extCharSpace = 16;

bool startsWithString(std::string line, std::string token) {
    return (line.substr(0, token.size()) == token);
}

}  // namespace

const char* getCellName(int nbDimension, int nbPoints) {
    if (nbDimension == 2) {
        if (nbPoints == 3) return triName;
        if (nbPoints == 4) return quadName;
    }
    if (nbDimension == 3) {
        if (nbPoints == 4) return tetraName;
        if (nbPoints == 5) return pyraName;
        if (nbPoints == 8) return hexaName;
    }

    return "UNDEFINED";
}

std::vector<int> findUniquePIDs(std::array<std::vector<int>, 4> array) {
    std::vector<int> vec;
    for (auto&& v : array) {
        vec.insert(vec.end(), v.begin(), v.end());
    }

    if (vec.size() == 0)
        vec = {defaultPID};
    else {
        std::sort(vec.begin(), vec.end());
        auto it = std::unique(vec.begin(), vec.end());
        vec.resize(std::distance(vec.begin(), it));
    }

    return vec;
}

D3D_status writeBDFFromCommonMeshData(CommonMeshData& mesh,
                            const boost::filesystem::path& path) {
    if (mesh.gridPoints.size() == 0) {
        std::cout << "No vertices to write\n";
        return D3D_status::FAIL;
    }

    int nMaterialId = 1;

    FILE* fp = NULL;
    if (path.string().c_str() == nullptr) return D3D_status::FAIL;
    fp = fopen(path.string().c_str(), "w");
    if (fp == NULL) return D3D_status::FAIL;

    // d3d header
    fprintf(fp,
            "TITLE = written by Divergent 3D system v1.0 (May 10 13:08:56"
            " 2018) \n ");
    fprintf(fp, "BEGIN BULK\n");

    fprintf(fp, "MAT1    %-8d%-8.1e        %-8f%-8.1e%-8.1e\n", nMaterialId,
            3e7, 0.33, 6.5e-6, 5.37e2);

    // if this is not empty, we need to make the connectivity refer to these ids
    if (mesh.gridIds.size() == 0) {
        std::cout << "Creating new indexing for the mesh grid points\n";
        for (auto ii = 0; ii < mesh.gridPoints.size(); ++ii) {
            mesh.gridIds.push_back(ii);
        }
    }

    auto pids = findUniquePIDs(mesh.cellPIDs);
    std::for_each(pids.begin(), pids.end(), [&](int pid) {
        fprintf(fp, "PSHELL  %-8d%-8d%-8f\n", pid, nMaterialId, 0.01 * pid);
    });

    // !! The bdf file format starts indexing at 1, so we need to add 1 to each
    // grid index
    for (auto dim = 0; dim < mesh.connectivity.size(); ++dim) {
        if (mesh.cellIds[dim].size() == 0) {
            for (auto ii = 0; ii < mesh.connectivity[dim].size(); ++ii) {
                mesh.cellIds[dim].push_back(ii + 1);
            }
        }

        if (mesh.connectivity[dim].size() > 0) {
            for (auto ii = 0; ii < mesh.connectivity[dim].size(); ++ii) {
                auto cellPID = mesh.cellPIDs[dim].size() > 0
                                   ? mesh.cellPIDs[dim][ii]
                                   : defaultPID;

                auto cell = mesh.connectivity[dim][ii];
                auto elemName = getCellName(dim, cell.size());
                fprintf(fp, "%s%-8d%-8d", elemName, mesh.cellIds[dim][ii],
                        cellPID);
                for (auto jj = 0; jj < cell.size(); ++jj) {
                    if (jj + 4 % 10 == 0) fprintf(fp, "\n        ");
                    fprintf(fp, "%-8d", cell[jj]);
                }
                fprintf(fp, "\n");
            }
        }
    }
    fprintf(fp, "PMASS   4       0.100000\n");

    for (auto ii = 0; ii < mesh.gridPoints.size(); ++ii) {
        auto coords = mesh.gridPoints[ii];
        fprintf(fp, "GRID*   %-32d%-16e%-16e*\n*       %-16e\n",
                mesh.gridIds[ii], coords[0], coords[1], coords[2]);
    }
    fprintf(fp, "ENDDATA\n");
    fclose(fp);
    return D3D_status::SUCCESS;
}

D3D_status readGridPoints(std::ifstream& femFile, std::string line,
                          CommonMeshData& mesh) {
    std::string _;
    int countPoints = 0;
    const int nCoords = 3;
    std::array<double, nCoords> coords;
    int gridId;
    bool continueReading = true;
    try {
        while (continueReading) {
            if (startsWithString(line, "GRID ")) {
                gridId = atoi(line.substr(8, standardCharSpace).c_str());
                coords[0] = atof(line.substr(24, standardCharSpace).c_str());
                coords[1] = atof(line.substr(32, standardCharSpace).c_str());
                coords[2] = atof(line.substr(40, standardCharSpace).c_str());
                mesh.gridPoints.push_back(coords);
                mesh.gridIds.push_back(gridId);
                countPoints++;
            } else if (startsWithString(line, "GRID*")) {
                gridId = atoi(line.substr(8, extCharSpace).c_str());
                coords[0] = atof(line.substr(40, extCharSpace).c_str());
                coords[1] = atof(line.substr(56, extCharSpace).c_str());
                getline(femFile, line);
                coords[2] = atof(line.substr(8, extCharSpace).c_str());
                mesh.gridPoints.push_back(coords);
                mesh.gridIds.push_back(gridId);
                countPoints++;
            } else if (startsWithString(line, "GRID")) {
                std::string field;
                std::istringstream iss(line);
                getline(iss, field, ',');  // GRID
                getline(iss, field, ',');  // point ID
                gridId = atoi(field.c_str());
                getline(iss, field, ',');  // other ID

                for (auto ii = 0; ii < nCoords; ++ii) {
                    getline(iss, field, ',');
                    coords[ii] = atof(field.c_str());
                }

                mesh.gridPoints.push_back(coords);
                mesh.gridIds.push_back(gridId);
                countPoints++;
            } else {
                break;
            }
            continueReading = getline(femFile, line).good();
        }
    } catch (...) {
        return D3D_status::CANNOT_READ_MESH;
    }

    return D3D_status::SUCCESS;
}

D3D_status readGridCells(
    std::ifstream& femFile, std::string line, d3d::CommonMeshData& mesh,
    std::map<int, int>& pointMap,
    std::array<std::vector<CommonMeshData::Cell>, 4>& tempConnectivity) {
    std::array<int, 4> cellNumbers;
    for (auto ii = 0; ii < cellNumbers.size(); ++ii) {
        cellNumbers[ii] = mesh.cellIds[ii].size();
    }

    bool continueReading = true;

    try {
        while (continueReading) {
            int charId = 8;
            bool ext;
            CellElement thisElemType;
            if (line[0] == '+') {
                continueReading = getline(femFile, line).good();
                continue;
            }
            if (line[0] != 'C') {
                break;
            }
            for (auto elemType : allCellElements) {
                if (startsWithString(line, elemType.name)) {
                    thisElemType = elemType;
                    ext = line[elemType.size] == '*';
                    break;
                }
            }

            if (thisElemType.defined) {
                int charSpace;
                charSpace = ext ? extCharSpace : standardCharSpace;

                auto cellId = atoi(line.substr(charId, charSpace).c_str());
                charId += charSpace;
                auto cellTag = atoi(line.substr(charId, charSpace).c_str());
                charId += charSpace;

                auto tempCellConnectivity =
                    CommonMeshData::Cell(thisElemType.numPoints);

                for (int ii = 0; ii < thisElemType.numPoints; ++ii) {
                    if (charId >= 72) {
                        getline(femFile,
                                line);  // skip to second line of coordinates
                        charId = 8;
                    }

                    auto gridId = atoi(line.substr(charId, charSpace).c_str());
                    tempCellConnectivity[ii] = gridId;
                    charId += charSpace;
                }

                tempConnectivity[thisElemType.dim].push_back(
                    tempCellConnectivity);

                mesh.cellIds[thisElemType.dim].push_back(cellId);
                mesh.cellPIDs[thisElemType.dim].push_back(cellTag);
                mesh.cellTypes[thisElemType.dim].push_back(thisElemType.type);

                ++cellNumbers[thisElemType.dim];
            } else {
                break;
            }

            continueReading = getline(femFile, line).good();
        }
    } catch (...) {
        return D3D_status::CANNOT_READ_MESH;
    }

    return D3D_status::SUCCESS;
}
D3D_status readBDFToCommonMeshData(const boost::filesystem::path& meshPath,
                         CommonMeshData& mesh) {
    int _;
    return readBDFToCommonMeshData(meshPath, mesh, _);
}

D3D_status readBDFToCommonMeshData(const boost::filesystem::path& meshPath,
                         CommonMeshData& mesh, int& designDomainPID) {
    std::ifstream femFile(meshPath.string());

    std::cout << "Read input file." << std::endl;

    std::map<int, int> pointMap;
    std::array<std::vector<CommonMeshData::Cell>, 4> tempConnectivity;

    auto ret_code = D3D_status::SUCCESS;

    if (femFile.is_open()) {
        std::string line;
        std::string _;

        while (getline(femFile, line)) {
            if (startsWithString(line, "DTPL")) {
                std::istringstream iss(line);
                iss >> _ >> _ >> _ >> designDomainPID;
                if (iss.fail()) {
                    iss.clear();
                    std::cout << "Design domain key could not be read."
                              << std::endl;
                    designDomainPID = -1;
                }
            }
            if (startsWithString(line, "GRID")) {
                ret_code = readGridPoints(femFile, line, mesh);
                if (ret_code != D3D_status::SUCCESS) {
                    return ret_code;
                }
                for (auto ii = 0; ii < mesh.gridIds.size(); ++ii) {
                    pointMap[mesh.gridIds[ii]] = ii;
                }
            }

            if (std::any_of(allCellElements.begin(), allCellElements.end(),
                            [line](CellElement e) {
                                return startsWithString(line, e.name);
                            })) {
                // we have no guarantee that the grid Ids have been read before.
                // Thats why we create a temp connectivity with the grid Ids,
                // and we create the real connectivity after having read
                // everything with the point map and the temp connectivity
                ret_code = readGridCells(femFile, line, mesh, pointMap,
                                         tempConnectivity);
                if (ret_code != D3D_status::SUCCESS) {
                    return ret_code;
                }
            }
        }

        auto numCells = 0;
        std::for_each(
            tempConnectivity.begin(), tempConnectivity.end(),
            [&](std::vector<std::vector<int>> vec) { numCells += vec.size(); });

        std::cout << "Finished reading " << mesh.gridPoints.size()
                  << " points and " << numCells << " cells." << std::endl;
    }

    // the tempConnectivity data structure refers to the vertices by their fem
    // grid ids. We use the point map to instead make them refer to their ids in
    // the new vectors

    for (auto ndim = 0; ndim < tempConnectivity.size(); ++ndim) {
        mesh.connectivity[ndim] =
            std::vector<std::vector<int>>(tempConnectivity[ndim].size());
        for (auto ii = 0; ii < tempConnectivity[ndim].size(); ++ii) {
            auto newCell = std::vector<int>(tempConnectivity[ndim][ii].size());
            for (auto jj = 0; jj < tempConnectivity[ndim][ii].size(); ++jj) {
                newCell[jj] = pointMap[tempConnectivity[ndim][ii][jj]];
            }
            mesh.connectivity[ndim][ii] = newCell;
        }
    }

    return D3D_status::SUCCESS;
}

}  // namespace io

}  // namespace d3d
