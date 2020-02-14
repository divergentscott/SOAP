
#include "meshData.h"

bool isIntInVector(int x, std::vector<int> vec) {
    auto it = std::find(vec.begin(), vec.end(), x);
    return (it != vec.end()) ? true : false;
}

int extractMeshDataByPID(const CommonMeshData &mesh_in, std::vector<int> pids,
                                CommonMeshData &mesh_out) {
    // mesh_out.gridIds = mesh_in.gridIds;
    mesh_out.gridPoints = mesh_in.gridPoints;

    bool pid_found = false;

    try {
        for (auto dim = 0; dim < mesh_in.connectivity.size(); ++dim) {
            for (auto ii = 0; ii < mesh_in.connectivity[dim].size(); ++ii) {
                if (mesh_in.cellPIDs[dim].size() > ii)
                    if (isIntInVector(mesh_in.cellPIDs[dim][ii], pids)) {
                        pid_found = true;
                        if (mesh_in.cellIds[dim].size() > ii)
                            mesh_out.cellIds[dim].push_back(
                                mesh_in.cellIds[dim][ii]);
                        mesh_out.cellPIDs[dim].push_back(
                            mesh_in.cellPIDs[dim][ii]);
                        mesh_out.connectivity[dim].push_back(
                            mesh_in.connectivity[dim][ii]);
                    }
            }
        }
    } catch (...) {
        return 0;
    }

    if (!pid_found) return 2;

    return 1;
}

int overwriteMeshDataPIDs(CommonMeshData &mesh_to_modify,
                                 const CommonMeshData &source_mesh,
                                 const std::vector<int> &pids) {
    bool pid_found = false;

    try {
        auto new_cell_id = 0;
        for (auto cellIds : mesh_to_modify.cellIds) {
            if (cellIds.size() > 0) {
                auto maxCellId =
                    *(std::max_element(cellIds.begin(), cellIds.end()));
                new_cell_id = std::max(0, maxCellId + 1);
            }
        }
        auto it_to_max_grid_id = std::max_element(
            mesh_to_modify.gridIds.begin(), mesh_to_modify.gridIds.end());
        auto new_grid_id = *(it_to_max_grid_id) + 1;

        CommonMeshData new_mesh;
        new_mesh.gridIds = mesh_to_modify.gridIds;
        new_mesh.gridPoints = mesh_to_modify.gridPoints;

        // copy grid points from the source mesh and make a mapping between the
        // new grid ids and the original ones
        std::vector<int> new_grid_ids;
        for (auto ii = 0; ii < source_mesh.gridPoints.size(); ++ii) {
            new_mesh.gridPoints.push_back(source_mesh.gridPoints[ii]);
            new_mesh.gridIds.push_back(new_grid_id);
            new_grid_ids.push_back(new_grid_id);
            ++new_grid_id;
        }

        // copy all cells from mesh_to_modify that have a different PID than the
        // input PIDs
        for (auto dim = 0; dim < mesh_to_modify.connectivity.size(); ++dim) {
            for (auto ii = 0; ii < mesh_to_modify.connectivity[dim].size();
                 ++ii) {
                if (mesh_to_modify.cellPIDs[dim].size() > ii)
                    if (!isIntInVector(mesh_to_modify.cellPIDs[dim][ii],
                                       pids)) {
                        if (mesh_to_modify.cellIds[dim].size() > ii)
                            new_mesh.cellIds[dim].push_back(
                                mesh_to_modify.cellIds[dim][ii]);
                        new_mesh.cellPIDs[dim].push_back(
                            mesh_to_modify.cellPIDs[dim][ii]);

                        std::vector<int> newCell;
                        for (auto nodeId :
                             mesh_to_modify.connectivity[dim][ii]) {
                            newCell.push_back(mesh_to_modify.gridIds[nodeId]);
                        }
                        new_mesh.connectivity[dim].push_back(newCell);

                        // new_mesh.connectivity[dim].push_back(
                        //    mesh_to_modify.connectivity[dim][ii]);
                    }
            }
        }

        // copy all cells from source_mesh that have a PID included in the input
        // PIDs. Use the new grid ids!

        for (auto dim = 0; dim < source_mesh.connectivity.size(); ++dim) {
            for (auto ii = 0; ii < source_mesh.connectivity[dim].size(); ++ii) {
                if (source_mesh.cellPIDs[dim].size() > ii)
                    if (isIntInVector(source_mesh.cellPIDs[dim][ii], pids)) {
                        pid_found = true;
                        if (source_mesh.cellIds[dim].size() > ii) {
                            new_mesh.cellIds[dim].push_back(new_cell_id);
                            new_cell_id++;
                        }
                        new_mesh.cellPIDs[dim].push_back(
                            source_mesh.cellPIDs[dim][ii]);

                        std::vector<int> newCell;
                        for (auto nodeId : source_mesh.connectivity[dim][ii]) {
                            newCell.push_back(new_grid_ids[nodeId]);
                        }
                        new_mesh.connectivity[dim].push_back(newCell);
                    }
            }
        }

        mesh_to_modify = new_mesh;

    } catch (...) {
        return 0;
    }

    if (!pid_found) return 2;

    return 1;
}

