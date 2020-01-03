//STL inclusions
#include <array>
#include <iostream>
#include <string>

//boost inclusions
#include "boost/filesystem.hpp"

//VTK inclusions
#include <vtkCenterOfMass.h>
#include <vtkClipPolyData.h>
#include <vtkContourTriangulator.h>
#include <vtkCutter.h>
#include <vtkOBBTree.h>
#include <vtkPlane.h>
#include <vtkPoints2D.h>
#include <vtkPolyData.h>
#include <vtkPolyDataNormals.h>
#include <vtkSmartPointer.h>
#include <vtkSTLReader.h>
#include <vtkStripper.h>
#include <vtkTriangle.h>
#include <vtkUnstructuredGrid.h>
#include <vtkXMLPolyDataWriter.h>

class Seam {

private:
	vtkSmartPointer<vtkPlane> plane_;
	vtkSmartPointer<vtkPolyData> curves_;
	vtkSmartPointer<vtkPolyData> disks_;
	bool is_disks_computed_ = false;
	double area_;
	bool is_area_computed_ = false;
	double perimeter_;
	bool is_perimeter_computed_ = false;
	bool is_mesh_sides_retained_ = false;
	vtkSmartPointer<vtkPolyData> above_;
	std::array<double, 3> above_box_size_{0.0,0.0,0.0 };
	vtkSmartPointer<vtkPolyData> below_;
	std::array<double, 3> below_box_size_{ 0.0,0.0,0.0 };

public:
	Seam(vtkSmartPointer<vtkPolyData> surface, vtkSmartPointer<vtkPlane> plane, bool retain_mesh_sides=false) {
		is_mesh_sides_retained_ = retain_mesh_sides;
		auto cutter = vtkSmartPointer<vtkCutter>::New();
		plane_ = plane;
		cutter->SetCutFunction(plane_);
		cutter->SetInputData(surface);
		cutter->Update();
		// Compute plane intersect surface
		curves_ = cutter->GetOutput();
		// Optimization potential: I'm running this extraction several times with different values. Very redundent.
		auto clipper = vtkSmartPointer<vtkClipPolyData>::New();
		clipper->SetClipFunction(plane_);
		clipper->SetInputData(surface);
		clipper->GenerateClippedOutputOn();
		clipper->Update();

		// Compute surface pieces cut by plane
		auto above = clipper->GetOutput();
		std::cout << "Above has:" << above->GetNumberOfPoints()<<"\n";
		double _corner[3], maxdir[3], middir[3], mindir[3], _sizes[3];
		auto obber = vtkSmartPointer<vtkOBBTree>::New();
		obber->ComputeOBB(above, _corner, maxdir, middir, mindir, _sizes);
		above_box_size_[0] = vtkMath::Norm(maxdir);
		above_box_size_[1] = vtkMath::Norm(middir);
		above_box_size_[2] = vtkMath::Norm(mindir);


		auto below = clipper->GetClippedOutput();
		std::cout << "Below has:" << below->GetNumberOfPoints() << "\n";
		obber->ComputeOBB(below, _corner, maxdir, middir, mindir, _sizes);
		below_box_size_[0] = vtkMath::Norm(maxdir);
		below_box_size_[1] = vtkMath::Norm(middir);
		below_box_size_[2] = vtkMath::Norm(mindir);

		if (is_mesh_sides_retained_) {
			above_ = above;
			below_ = below;
		}

	};

	vtkSmartPointer<vtkPolyData> get_curves() {
		return curves_;
	};

	std::array<double, 3> get_above_box_size() {
		return above_box_size_;
	}

	std::array<double, 3> get_below_box_size() {
		return below_box_size_;
	}

	void compute_perimeter() {
		double perimeter = 0;
		for (int edge_i = 0; edge_i < curves_->GetNumberOfCells(); edge_i++) {
			auto edge = curves_->GetCell(edge_i);
			double pt_coords[2][3];
			for (int foo = 0; foo < 2; foo++) edge->GetPoints()->GetPoint(foo, pt_coords[foo]);
			perimeter += sqrt(vtkMath::Distance2BetweenPoints(pt_coords[0], pt_coords[1]));
		}
		perimeter_ = perimeter;
		is_perimeter_computed_ = true;
	}

	double get_perimeter() {
		if (!is_perimeter_computed_) {
			compute_perimeter();
		}
		return perimeter_;
	}

	void compute_disks() {
		auto filler = vtkSmartPointer<vtkContourTriangulator>::New();
		filler->SetInputData(curves_);
		filler->Update();
		disks_ = filler->GetOutput();
		is_disks_computed_ = true;
	}

	vtkSmartPointer<vtkPolyData> get_disks() {
		if (!is_disks_computed_) {
			compute_disks();
		}
		return disks_;
	}

	void compute_area() {
		if (!is_disks_computed_) {
			compute_disks();
		}
		double area = 0;
		for (int tri_i = 0; tri_i < disks_->GetNumberOfCells(); tri_i++) {
			auto tri = disks_->GetCell(tri_i);
			double pt_coords[3][3];
			for (int foo = 0; foo < 3; foo++) tri->GetPoints()->GetPoint(foo, pt_coords[foo]);
			area += vtkTriangle::TriangleArea(pt_coords[0], pt_coords[1], pt_coords[2]);
		}
		area_ = area;
		is_area_computed_ = true;
	}

	double get_area() {
		if (!is_area_computed_) {
			compute_area();
		}
		return area_;
	}

};

class CrossSection {
private:
	vtkSmartPointer<vtkPlane> plane_;
	//vtkSmartPointer<vtkPolyData> polygons_ = vtkSmartPointer<vtkPolyData>::New();
	double planarx_[3] { 0.0, 0.0, 0.0 };
	double planary_[3]{ 0.0, 0.0, 0.0 };
	double cut_origin_[3]{ 0.0, 0.0, 0.0 };

	CrossSection(vtkSmartPointer<vtkPlane> plane, vtkSmartPointer<vtkPolyData> disks) {
		plane_ = plane;
		double planenormal[3];
		plane_->GetNormal(planenormal);
		vtkMath::Normalize(planenormal);
		double planeorigin[3];
		plane_->GetOrigin(planeorigin);
		//Compute cut origin
		double origin_tmp0[3];
		vtkCenterOfMass::ComputeCenterOfMass(disks->GetPoints(), nullptr, origin_tmp0);
		plane_->ProjectPoint(origin_tmp0, cut_origin_);
		//The ppoints are the in-plane points: planar points.
		auto ppoints = vtkSmartPointer<vtkPoints2D>::New();
		//Project all points of the cross section polys into the plane with the {planarx_, planary_} basis
	}
};

class SeameyTodd {
	//SeameyTodd slices his patients into pieces and Mrs. Cuttit will tell us how they taste
private:
	vtkSmartPointer<vtkPolyData> patient_; 
public:
	Seam slice(vtkSmartPointer<vtkPlane> plane) {
		return Seam(patient_, plane);
	}
};

int main() {
	boost::filesystem::path infilepath{"C:\\Users\\sscott\\Pictures\\trex_connected.stl"};
	auto reader = vtkSmartPointer<vtkSTLReader>::New();
	reader->SetFileName(infilepath.string().c_str());
	reader->Update();
	auto normaler = vtkSmartPointer<vtkPolyDataNormals>::New();
	normaler->SetInputConnection(reader->GetOutputPort());
	normaler->ComputeCellNormalsOn();
	normaler->ComputePointNormalsOn();
	normaler->Update();
	auto surface = normaler->GetOutput();



	std::cout << surface->GetNumberOfPoints();
	auto sharpplane = vtkSmartPointer<vtkPlane>::New();
	double anorigin[3]{ 100.0, 100.0, 10.0 };
	sharpplane->SetOrigin(anorigin);
	double adirection[3]{ 0.1235, 0.1235, 0.1235};
	sharpplane->SetNormal(adirection);

	auto seamly = Seam(surface, sharpplane);
	auto curves = seamly.get_curves();
	std::cout << curves->GetNumberOfPoints() << " " << curves->GetNumberOfCells() << " " << curves->GetNumberOfLines() << "\n";

	auto writer = vtkSmartPointer<vtkXMLPolyDataWriter>::New();
	writer->SetInputData(curves);
	writer->SetFileName("aquickcurve.vtp");
	writer->Write();

	auto writer2 = vtkSmartPointer<vtkXMLPolyDataWriter>::New();
	writer2->SetInputData(seamly.get_disks());
	writer2->SetFileName("aquickdisk.vtp");
	writer2->Write();

	std::cout << "Perimeter: " << seamly.get_perimeter() << "\n";
	std::cout << "Area: " << seamly.get_area() << "\n";

	std::cout << "Above size:";
	for (auto x : seamly.get_above_box_size()) std::cout << " " << x;
	std::cout<< "\n";

	std::cout << "Below size:";
	for (auto x : seamly.get_below_box_size()) std::cout << " " << x;
	std::cout << "\n";

};