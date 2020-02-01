
#include "meshData.h"

namespace d3d {

bool isIntInVector(int x, std::vector<int> vec) {
    auto it = std::find(vec.begin(), vec.end(), x);
    return (it != vec.end()) ? true : false;
}

D3D_status extractMeshDataByPID(const CommonMeshData &mesh_in, std::vector<int> pids,
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
        return D3D_status::FAIL;
    }

    if (!pid_found) return D3D_status::PID_NOT_FOUND;

    return D3D_status::SUCCESS;
}

D3D_status overwriteMeshDataPIDs(CommonMeshData &mesh_to_modify,
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
        return D3D_status::FAIL;
    }

    if (!pid_found) return D3D_status::PID_NOT_FOUND;

    return D3D_status::SUCCESS;
}


vtkSmartPointer<vtkPolyData> makePolydata(const CommonMeshData& mesh) {
	vtkSmartPointer<vtkPoints> points = vtkSmartPointer<vtkPoints>::New();
	auto polydata = vtkSmartPointer<vtkPolyData>::New();
	std::for_each(
		mesh.gridPoints.begin(), mesh.gridPoints.end(),
		[&](CommonMeshData::Point pt) { points->InsertNextPoint(pt.data()); });
	polydata->SetPoints(points);
	vtkSmartPointer<vtkIdTypeArray> connectivityArray =
		vtkSmartPointer<vtkIdTypeArray>::New();
	std::for_each(mesh.connectivity[2].begin(), mesh.connectivity[2].end(),
		[&](CommonMeshData::Cell cell) {
		connectivityArray->InsertNextValue(3);
		std::for_each(cell.begin(), cell.end(), [&](int id) {
			connectivityArray->InsertNextValue(id);
		});
	});
	vtkSmartPointer<vtkCellArray> cellArray =
		vtkSmartPointer<vtkCellArray>::New();
	cellArray->SetCells(mesh.connectivity[2].size(), connectivityArray);
	polydata->SetPolys(cellArray);
	return polydata;
}

CommonMeshData makeCommonMeshData(vtkSmartPointer<vtkPolyData> polydata) {
	CommonMeshData mesh;
	for (auto ii = 0; ii < polydata->GetNumberOfPoints(); ++ii) {
		CommonMeshData::Point point = { polydata->GetPoint(ii)[0],
										polydata->GetPoint(ii)[1],
										polydata->GetPoint(ii)[2] };
		mesh.gridPoints.push_back(point);
		mesh.gridIds.push_back(ii);
	}
	for (auto jj = 0; jj < polydata->GetNumberOfCells(); ++jj) {
		auto cellPointIds = vtkSmartPointer<vtkIdList>::New();
		polydata->GetCellPoints(jj, cellPointIds);
		CommonMeshData::Cell cell;
		for (auto kk = 0; kk < cellPointIds->GetNumberOfIds(); ++kk) {
			cell.push_back(cellPointIds->GetId(kk));
		}
		auto dim = polydata->GetCellType(jj) == VTK_QUAD ||
			polydata->GetCellType(jj) == VTK_TRIANGLE
			? 2
			: 3;
		mesh.connectivity[dim].push_back(cell);
		mesh.cellPIDs[dim].push_back(1);
		mesh.cellIds[dim].push_back(jj);
	}
	return mesh;
}
CommonMeshData makeCommonMeshData(vtkSmartPointer<vtkUnstructuredGrid> uGrid) {
	CommonMeshData mesh;
	mesh.gridPoints.resize(uGrid->GetNumberOfPoints());
	for (auto ii = 0; ii < uGrid->GetNumberOfPoints(); ++ii) {
		auto point = uGrid->GetPoint(ii);
		mesh.gridPoints[ii] = { point[0], point[1], point[2] };
	}
	mesh.connectivity[3].resize(uGrid->GetNumberOfCells());
	for (auto ii = 0; ii < uGrid->GetNumberOfCells(); ++ii) {
		auto cell = uGrid->GetCell(ii);
		mesh.connectivity[3][ii].resize(cell->GetNumberOfPoints());
		for (auto jj = 0; jj < cell->GetNumberOfPoints(); ++jj) {
			mesh.connectivity[3][ii][jj] = cell->GetPointId(jj);
		}
	}
	return mesh;
};
vtkSmartPointer<vtkUnstructuredGrid> makeUGrid(
	const d3d::CommonMeshData& mesh) {
	vtkSmartPointer<vtkUnstructuredGrid> uGrid =
		vtkSmartPointer<vtkUnstructuredGrid>::New();
	vtkSmartPointer<vtkPoints> points = vtkSmartPointer<vtkPoints>::New();
	std::for_each(mesh.gridPoints.begin(), mesh.gridPoints.end(),
		[&](CommonMeshData::Point coords) {
		points->InsertNextPoint(coords.data());
	});
	uGrid->SetPoints(points);
	vtkSmartPointer<vtkIdTypeArray> connectivityArray =
		vtkSmartPointer<vtkIdTypeArray>::New();
	std::vector<int> cellTypes;
	for (int ii = 0; ii < mesh.connectivity[3].size(); ii++) {
		auto cell = mesh.connectivity[3][ii];
		auto cellSize = cell.size();
		connectivityArray->InsertNextValue(cellSize);
		std::for_each(cell.begin(), cell.end(), [&](int pointId) {
			connectivityArray->InsertNextValue(pointId);
		});
		if (cellSize == 4)
			cellTypes.push_back(VTK_TETRA);
		else if (cellSize == 5)
			cellTypes.push_back(VTK_PYRAMID);
		else
			cellTypes.push_back(VTK_HEXAHEDRON);
	}
	vtkSmartPointer<vtkCellArray> cellArray =
		vtkSmartPointer<vtkCellArray>::New();
	cellArray->SetCells(cellTypes.size(), connectivityArray);
	uGrid->SetCells(cellTypes.data(), cellArray);
	return uGrid;
}


}  // namespace d3d
