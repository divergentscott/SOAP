#include <iostream>

#include "bdfIO.h"
#include "seamDesigner.h"

#include <vtkAppendPolyData.h>
#include <vtkLine.h>
#include <vtkSTLReader.h>
#include <vtkTransform.h>
#include <vtkTransformPolyDataFilter.h>
#include <vtkTriangle.h>
#include <vtkXMLPolyDataWriter.h>


// Tao is 2pi
const double CONSTANT_TAO = 6.28318530717958647692528676655900576839433879875021164194;


void write_polydata(const vtkSmartPointer<vtkPolyData> x, const std::string filename) {
	//Write
	auto writer3 = vtkSmartPointer<vtkXMLPolyDataWriter>::New();
	writer3->SetInputData(x);
	writer3->SetFileName(filename.c_str());
	writer3->Write();
}


vtkSmartPointer<vtkPolyData> generate_an_ngon(const int nn) {
	auto ngon_polydata = vtkSmartPointer<vtkPolyData>::New();
	auto ngon_pts = vtkSmartPointer<vtkPoints>::New();
	std::vector<std::vector<int>> cells;
	for (int foo = 0; foo < nn; foo++) {
		double ptcoor[3]{ cos(CONSTANT_TAO / nn * foo), sin(CONSTANT_TAO / nn * foo), 0.0 };
		ngon_pts->InsertPoint(foo, ptcoor);
		//std::vector<int> cell
		cells.push_back({ foo, (foo + 1)%nn });
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

//void example_check_reader() {
//	d3d::CommonMeshData meshin;
//	boost::filesystem::path filein = "C:\\Users\\sscott\\Pictures\\unitsphere_meshlab.bdf";
//	auto err = d3d::io::readBDFToCommonMeshData(filein, meshin);
//	d3d::soap::Designer::Parameters paramers;
//	paramers.plane_origin = { 462.15894026125324, 119.3025104782025, -100.30062516465588 };
//	paramers.plane_normal = { -0.1398971605146768, 0.9816945964030201, -0.12924590466642394 };
//	paramers.tongue_direction = { 0.0, 1.0, 0.0 };
//	paramers.groove_outer = 0.01;
//	paramers.groove_inner = -0.01;
//	paramers.gap_radial = 0.005;
//	paramers.gap_depth = 0.005;
//	paramers.tongue_depth = 0.05;
//	paramers.trim_depth = 0.025;
//	paramers.use_tongue = true;
//	d3d::soap::Designer dessy = d3d::soap::Designer(meshin, paramers);
//}

void example_check_pentagon() {
	auto penta = generate_an_ngon(5);
	write_polydata(penta, "started_a_penta.vtp");
	auto curve = d3d::planar::CurveCollection(penta);
	curve.write_to_vtp("perhaps_a_pentagon.vtp");
}

void example_shape_parade() {
	auto appender = vtkSmartPointer<vtkAppendPolyData>::New();

	auto s3 = generate_an_ngon(3);
	appender->AddInputData(s3);

	auto a4 = generate_an_ngon(3);
	auto transform4 = vtkSmartPointer<vtkTransform>::New();
	auto transformer4 = vtkSmartPointer<vtkTransformPolyDataFilter>::New();
	transform4->Scale(2.0, 2.0, 2.0);
	transformer4->SetTransform(transform4);
	transformer4->SetInputData(a4);
	transformer4->Update();
	appender->AddInputConnection(transformer4->GetOutputPort());
	
	auto a5 = generate_an_ngon(5);
	auto transform5 = vtkSmartPointer<vtkTransform>::New();
	auto transformer5 = vtkSmartPointer<vtkTransformPolyDataFilter>::New();
	transform5->Translate(3.0, 0.0, 0.0);
	transformer5->SetTransform(transform5);
	transformer5->SetInputData(a5);
	transformer5->Update();
	appender->AddInputConnection(transformer5->GetOutputPort());

	auto a6 = generate_an_ngon(6);
	auto transform6 = vtkSmartPointer<vtkTransform>::New();
	auto transformer6 = vtkSmartPointer<vtkTransformPolyDataFilter>::New();
	transform6->Translate(0.0, 3.0, 0.0);
	transformer6->SetTransform(transform6);
	transformer6->SetInputData(a6);
	transformer6->Update();
	appender->AddInputConnection(transformer6->GetOutputPort());
	
	appender->Update();
	auto curve = d3d::planar::CurveCollection(appender->GetOutput());
	curve.write_to_vtp("shape_parade.vtp");
}


vtkSmartPointer<vtkPolyData> scaledown(vtkSmartPointer<vtkPolyData> x, const double scale) {
	auto transform4 = vtkSmartPointer<vtkTransform>::New();
	auto transformer4 = vtkSmartPointer<vtkTransformPolyDataFilter>::New();
	transform4->Scale(scale, scale, scale);
	transformer4->SetTransform(transform4);
	transformer4->SetInputData(x);
	transformer4->Update();
	return transformer4->GetOutput();
};

void example_lamination(int iterates) {
	auto appender = vtkSmartPointer<vtkAppendPolyData>::New();
	auto s3 = generate_an_ngon(5);
	appender->AddInputData(s3);
	std::vector<vtkSmartPointer<vtkPolyData>> polydatas;
	for (int foo = 1; foo < iterates; foo++) {
		double scale = 1 + foo * 0.0001;
		polydatas.push_back(scaledown(s3,scale));
	}
	for (auto x : polydatas) { appender->AddInputData(x); };
	appender->Update();
	auto lamin = appender->GetOutput();
	auto curve = d3d::planar::CurveCollection(lamin);
	curve.write_to_vtp("lamination.vtp");
}

void example_offset() {
	auto s3 = generate_an_ngon(5);
	auto curve = d3d::planar::CurveCollection(s3);
	auto off = curve.distance_field(-0.1, 0.4, 23);
	write_polydata(off, "pentoff.vtp");
}

void example_cs() {
	auto reader = vtkSmartPointer<vtkSTLReader>::New();
	reader->SetFileName("C:\\Users\\sscott\\Pictures\\unitsphere_meshlab.stl");
	reader->Update();
	auto mesh = reader->GetOutput();
	d3d::soap::point::r3 normal{ 0.0, 0.0, 1.0 };
	d3d::soap::point::r3 origin{ 0.0, 0.0, 0.0 };
	d3d::soap::CrossSectioner cs(mesh, origin, normal);
	write_polydata(cs.get_cross_section(), "stillacircle1.vtp");
	write_polydata(cs.get_planed_cross_section(),"stillacircle2.vtp");
}

//void example_sphere_seam() {
//	auto reader = vtkSmartPointer<vtkSTLReader>::New();
//	reader->SetFileName("C:\\Users\\sscott\\Pictures\\unitsphere_meshlab.stl");
//	reader->Update();
//	auto mesh = reader->GetOutput();
//	d3d::soap::SeamDesigner::Parameters params {
//		{0.02}, //groove_outer
//		{-0.02}, //groove_inner
//		{0.01}, //gap_radial
//		{0.01}, //gap_depth
//		{0.06}, //tongue_depth
//		{0.02}, //trim_depth
//		{{ 0.0, 0.0, 0.0 }}, //plane_origin
//		{{ 0.0, 0.0, 1.0 }}, //plane_normal
//		{false} //use_tongue
//		//{{ 0.0, 0.0, 0.0 }}, //tongue_direction
//	};
//	auto seam_designer = d3d::soap::SeamDesigner(mesh, params);
//	write_polydata(seam_designer.get_seam(), "orb_sliced.vtp");
//}


int main() {
	example_cs();
};