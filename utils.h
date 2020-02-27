#ifndef MESHDATALIB_UTILS_H
#define MESHDATALIB_UTILS_H

#include <vtkPolyData.h>
#include <vtkSmartPointer.h>
#include <vtkUnstructuredGrid.h>

#include "meshData.h"

vtkSmartPointer<vtkPolyData> makePolydata(const CommonMeshData& mesh);
CommonMeshData makeCommonMeshData(vtkSmartPointer<vtkPolyData> polydata);

vtkSmartPointer<vtkUnstructuredGrid> makeUGrid(const CommonMeshData& mesh);
CommonMeshData makeCommonMeshData(vtkSmartPointer<vtkUnstructuredGrid> uGrid);

#endif