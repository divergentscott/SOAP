#include <iostream>

#include "bdfIO.h"
#include "seamDesigner.h"
#include <vtkLine.h>
#include <vtkTriangle.h>


// Tao is 2pi
const double CONSTANT_TAO = 6.28318530717958647692528676655900576839433879875021164194;

vtkSmartPointer<vtkPolyData> generate_an_ngon(const int nn) {
	auto ngon_polydata = vtkSmartPointer<vtkPolyData>::New();
	auto ngon_pts = vtkSmartPointer<vtkPoints>::New();
	std::vector<std::vector<int>> cells;
	for (int foo = 0; foo < nn; foo++) {
		double ptcoor[3]{ cos(CONSTANT_TAO / nn * foo), sin(CONSTANT_TAO / nn * foo), 0.0 };
		ngon_pts->InsertPoint(foo, ptcoor);
		//std::vector<int> cell
		if ((foo > 0)&(foo < (nn - 1))) cells.push_back({ foo, foo + 1});
	}
	auto ngon_cells = vtkSmartPointer<vtkCellArray>::New();
	// Copy the triangles to the polydata
	for (auto cell : cells) {
		auto lien = vtkSmartPointer<vtkLine>::New();
		for (int foo = 0; foo < 2; foo++) lien->GetPointIds()->SetId(foo, cell[foo]);
		ngon_cells->InsertNextCell(lien);
	}
	ngon_polydata->SetPoints(ngon_pts);
	ngon_polydata->SetLines(ngon_cells);
	return ngon_polydata;
}

void example_check_reader() {
	d3d::CommonMeshData meshin;
	boost::filesystem::path filein = "C:\\Users\\sscott\\Pictures\\unitsphere_meshlab.bdf";
	auto err = d3d::io::readBDFToCommonMeshData(filein, meshin);
	d3d::soap::Designer::Parameters paramers;
	paramers.plane_origin = { 462.15894026125324, 119.3025104782025, -100.30062516465588 };
	paramers.plane_normal = { -0.1398971605146768, 0.9816945964030201, -0.12924590466642394 };
	paramers.tongue_direction = { 0.0, 1.0, 0.0 };
	paramers.groove_outer = 0.01;
	paramers.groove_inner = -0.01;
	paramers.gap_radial = 0.005;
	paramers.gap_depth = 0.005;
	paramers.tongue_depth = 0.05;
	paramers.trim_depth = 0.025;
	paramers.use_tongue = true;
	d3d::soap::Designer dessy = d3d::soap::Designer(meshin, paramers);
}

int main() {
	std::cout << "Saluton Mundo!\n";
	auto penta =  generate_an_ngon(5);
	auto curve = d3d::planar::CurveCollection(penta);
};