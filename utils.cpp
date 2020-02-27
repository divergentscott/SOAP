#include "utils.h"

#include <vtkCellArray.h>
#include <vtkCellType.h>
#include <vtkIdTypeArray.h>
#include <vtkPoints.h>
#include <vtkPolyData.h>
#include <vtkSmartPointer.h>


vtkSmartPointer<vtkPolyData> makePolydata(const CommonMeshData& mesh) {
    vtkSmartPointer<vtkPoints> points = vtkSmartPointer<vtkPoints>::New();
    auto polydata = vtkSmartPointer<vtkPolyData>::New();

    std::for_each(
        mesh.gridPoints.begin(), mesh.gridPoints.end(),
        [&](CommonMeshData::Point pt) { points->InsertNextPoint(pt.data()); });
    polydata->SetPoints(points);

    vtkSmartPointer<vtkIdTypeArray> connectivityArray =
        vtkSmartPointer<vtkIdTypeArray>::New();

    for (auto cell : mesh.connectivity[2]) {
        connectivityArray->InsertNextValue(3);
        for (auto id : cell) {
            connectivityArray->InsertNextValue(id);
        }
    }

    vtkSmartPointer<vtkCellArray> cellArray =
        vtkSmartPointer<vtkCellArray>::New();
    cellArray->SetCells(mesh.connectivity[2].size(), connectivityArray);
    polydata->SetPolys(cellArray);

    return polydata;
}

CommonMeshData makeCommonMeshData(vtkSmartPointer<vtkPolyData> polydata) {
    CommonMeshData mesh;
    for (auto ii = 0; ii < polydata->GetNumberOfPoints(); ++ii) {
        CommonMeshData::Point point = {polydata->GetPoint(ii)[0],
                                       polydata->GetPoint(ii)[1],
                                       polydata->GetPoint(ii)[2]};
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
        mesh.gridPoints[ii] = {point[0], point[1], point[2]};
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
    const CommonMeshData& mesh) {
    vtkSmartPointer<vtkUnstructuredGrid> uGrid =
        vtkSmartPointer<vtkUnstructuredGrid>::New();
    vtkSmartPointer<vtkPoints> points = vtkSmartPointer<vtkPoints>::New();

	for (auto coords: mesh.gridPoints) points->InsertNextPoint(coords.data());

    uGrid->SetPoints(points);

    vtkSmartPointer<vtkIdTypeArray> connectivityArray =
        vtkSmartPointer<vtkIdTypeArray>::New();
    std::vector<int> cellTypes;

    for (int ii = 0; ii < mesh.connectivity[3].size(); ii++) {
        auto cell = mesh.connectivity[3][ii];
        auto cellSize = cell.size();
        connectivityArray->InsertNextValue(cellSize);

		for (auto pointId: cell) connectivityArray->InsertNextValue(pointId);

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

