#ifndef MESHDATALIB_MESHDATA_H
#define MESHDATALIB_MESHDATA_H

#include <array>
#include <boost/filesystem.hpp>
#include <vector>
#include "d3derr.h"

#include <vtkCellArray.h>
#include <vtkCellType.h>
#include <vtkIdTypeArray.h>
#include <vtkPoints.h>
#include <vtkPolyData.h>
#include <vtkSmartPointer.h>
#include <vtkUnstructuredGrid.h>

namespace d3d {

enum class CellType { CTRIA3, CQUAD4, CTETRA, CHEXA, CPYRA, UNDEFINED };

struct CellElement {
    std::string name;
    int size;
    int vtkType;
    CellType type = CellType::UNDEFINED;
    int numPoints;
    bool defined = false;
    int dim = -1;  // -1 means that the dimension is undefined

    CellElement(){};
    CellElement(std::string n, int vtkt, int nPts, int dim, CellType t)
        : name(n),
          size(n.size()),
          vtkType(vtkt),
          numPoints(nPts),
          dim(dim),
          type(t),
          defined(true){};
};

const std::vector<CellElement> allCellElements;

struct CommonMeshData {
    using Point = std::array<double, 3>;
    using Cell = std::vector<int>;

    std::vector<Point> gridPoints;
    std::vector<int> gridIds;

    // in the following structures, a cell will have its info in
    // connectivity[nDim][cellIndex][nodeIndex], cellPIDs[nDim][cellIndex],
    // cellIds[nDim][cellIndex]

    std::array<std::vector<Cell>, 4> connectivity;
    std::array<std::vector<int>, 4> cellIds;
    std::array<std::vector<int>, 4> cellPIDs;
    std::array<std::vector<CellType>, 4> cellTypes;
};

// mesh_out is a copy of mesh_in, but it only contains the cells with the
// specified pids
D3D_status extractMeshDataByPID(const CommonMeshData &mesh_in, std::vector<int> pids,
                                CommonMeshData &mesh_out);

// In the mesh to modify, all cells with the specified pids get overwritten.
// This mesh to modify keeps all other cells, and those cells and grid points
// keep their original ids. Typical use would be mesh_to_modify has design and
// non design space, source_mesh has the updated design space: after this
// function, mesh_to_modify keeps the non design space intact and has the
// modified design space
D3D_status overwriteMeshDataPIDs(CommonMeshData &mesh_to_modify,
                                 const CommonMeshData &source_mesh,
                                 const std::vector<int> &pids);


vtkSmartPointer<vtkPolyData> makePolydata(const CommonMeshData& mesh);
CommonMeshData makeCommonMeshData(vtkSmartPointer<vtkPolyData> polydata);
CommonMeshData makeCommonMeshData(vtkSmartPointer<vtkUnstructuredGrid> uGrid);
vtkSmartPointer<vtkUnstructuredGrid> makeUGrid(
	const d3d::CommonMeshData& mesh);

}  // namespace d3d

#endif
