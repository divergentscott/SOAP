//STL inclusions
#include <array>
#include <cstdlib>
#include <fstream>
#include <iostream>
#include <set>
#include <string>


#define _USE_MATH_DEFINES
#include <math.h>

//boost inclusions
#include "boost/filesystem.hpp"

//VTK inclusions
#include <vtkCenterOfMass.h>
#include <vtkCellArray.h>
#include <vtkClipPolyData.h>
#include <vtkContourTriangulator.h>
#include <vtkCutter.h>
#include <vtkDoubleArray.h>
#include <vtkFeatureEdges.h>
#include <vtkLine.h>
#include <vtkOBBTree.h>
#include <vtkPlane.h>
#include <vtkPoints2D.h>
#include <vtkPointData.h>
#include <vtkPolyData.h>
#include <vtkPolyDataNormals.h>
#include <vtkSmartPointer.h>
#include <vtkSTLReader.h>
#include <vtkStripper.h>
#include <vtkTransform.h>
#include <vtkTransformPolyDataFilter.h>
#include <vtkTriangle.h>
#include <vtkUnstructuredGrid.h>
#include <vtkVertex.h>
#include <vtkXMLPolyDataWriter.h>

using ppoint = std::array<double, 2>;

ppoint operator+(ppoint a, ppoint b) {
	ppoint c{ a[0] + b[0], a[1] + b[1] };
	return c;
}

ppoint operator-(ppoint a, ppoint b) {
	ppoint c{ a[0] - b[0], a[1] - b[1] };
	return c;
}

ppoint operator*(double a, ppoint b) {
	ppoint c{ a*b[0], a*b[1] };
	return c;
}


class PlanarNormalFilter {
private:
	vtkSmartPointer<vtkPolyData> polys_ = vtkSmartPointer<vtkPolyData>::New();
	vtkSmartPointer<vtkPolyData> boundary_ = vtkSmartPointer<vtkPolyData>::New();
	vtkSmartPointer<vtkDoubleArray> planar_normals_ = vtkSmartPointer<vtkDoubleArray>::New();
	double min_boundary_edge_length_{ -1.0 };
	const double repsilon_ = std::sqrt(std::numeric_limits<double>::epsilon());

public:
	PlanarNormalFilter(vtkSmartPointer<vtkPolyData> polys) {
		std::cout << "sqrt epsilon is: " << repsilon_ << "\n";
		polys_ = polys;
		auto edger = vtkSmartPointer<vtkFeatureEdges>::New();
		edger->BoundaryEdgesOn();
		edger->SetInputData(polys);
		edger->Update();
		auto boundary_edges = edger->GetOutput();
		boundary_ = boundary_edges;
		boundary_->BuildLinks();
		polys_->BuildLinks();
		//Get the minimal edge length. Used for normal orientation checking.
		for (int edge_i = 0; edge_i < boundary_->GetNumberOfCells(); edge_i++) {
			auto pts_in_i = boundary_->GetCell(edge_i)->GetPointIds();
			double pt0[3], pt1[3];
			boundary_->GetPoint(pts_in_i->GetId(0), pt0);
			boundary_->GetPoint(pts_in_i->GetId(1), pt1);
			double length_i = vtkMath::Distance2BetweenPoints(pt0, pt1);
			if (edge_i == 0) min_boundary_edge_length_ = length_i;
			else
			{
				min_boundary_edge_length_ = std::min(min_boundary_edge_length_, length_i);
			}
		}
	};
	
	std::vector<int> get_adj_boundary_points(int point_id) {
		auto ngb_cells = vtkSmartPointer<vtkIdList>::New();
		boundary_->GetPointCells(point_id, ngb_cells);
		std::set<int> ngb_points;
		for (int cell_i = 0; cell_i < ngb_cells->GetNumberOfIds(); cell_i++) {
			int cell = ngb_cells->GetId(cell_i);
			auto pts_in_i = boundary_->GetCell(cell)->GetPointIds();
			for (int pt_j = 0; pt_j < pts_in_i->GetNumberOfIds(); pt_j++) {
				int pt = pts_in_i->GetId(pt_j);
				if (pt != point_id) ngb_points.insert(pt);
			}
		}
		std::vector<int> ngb_points_vec(ngb_points.begin(), ngb_points.end());
		return ngb_points_vec;
	}

	void printit(double *data) {
		std::cout << data[0] << " " << data[1] << " " << data[2] <<  "\n";
	}

	vtkSmartPointer<vtkDoubleArray> get_raw_normals() {
		// Really the 2nd discrete path derivative of the bounding curve.
		// ~Normals except not normalized, and no consistency checking.
		auto raw_normals = vtkSmartPointer<vtkDoubleArray>::New();
		raw_normals->SetNumberOfComponents(2);
		raw_normals->SetNumberOfTuples(boundary_->GetNumberOfPoints());
		raw_normals->SetName("Raw Normals");
		for (int pt_i = 0; pt_i < boundary_->GetNumberOfPoints(); pt_i++) {
			//Get coordinates of the adjacent boundary points.
			std::vector<int> pts_adj = get_adj_boundary_points(pt_i);
			if (pts_adj.size() != 2) std::cout << "BOUNDARYISCRAZY!! @ " << pt_i << "\n";
			double pa[3], pb[3], pc[3], normal_i[3];
			boundary_->GetPoint(pts_adj[0], pa);
			boundary_->GetPoint(pts_adj[1], pc);
			boundary_->GetPoint(pt_i, pb);
			// For points a ~ b ~ c define the normal at b as proportional to a+c-2b
			double areaabc = (pc[0] - pb[0])*(pa[1] - pb[1]) - (pc[1] - pb[1])*(pa[0] - pb[0]);
			if (std::abs(areaabc) > repsilon_) {
				//Non-colinear abc
				for (int foo = 0; foo < 3; foo++) {
					normal_i[foo] = pa[foo] + pc[foo] - 2 * pb[foo];
				}
				vtkMath::Normalize2D(normal_i);
			} else {
				//Colinear
				normal_i[0] = pa[1] - pc[1];
				normal_i[1] = pc[0] - pa[0];
				normal_i[2] = 0;
				vtkMath::Normalize2D(normal_i);
			}
			raw_normals->SetTuple2(pt_i, normal_i[0], normal_i[1]);
		}
		return raw_normals;
	}

	std::array<double,2> get_tangent(int boundary_point_id) {
		std::vector<int> pts_adj = get_adj_boundary_points(boundary_point_id);
		double pa[3]{ 0.0,0.0,0.0 }, pc[3]{ 0.0,0.0,0.0 };
		boundary_->GetPoint(pts_adj[0], pa);
		boundary_->GetPoint(pts_adj[1], pc);
		// For points a ~ b ~ c define the tangent at b as proportional to c-a
		std::array<double, 2> tan;
		tan[0] = pc[0] - pa[0];
		tan[1] = pc[1] - pa[1];
		return tan;
	}

	std::array<double, 2> ray_segment_intersect(std::array<double, 2>orig, std::array<double, 2>r, std::array<double, 2>a, std::array<double, 2>b) {
		// Compute ray, line segment intersection.
		// Ray has origin orig and direction r. That is: { orig + t * r | t>0 } 
		// Line segment endpoints a b. That is: { a * s + b * (1-s) | s in [0,1] }
		// Computes the parameter t and s at intersection.
		// If the ray and segment are colinear, there may be infinitely many solutions.
		// In that case return the s = 0 solution, so t is such that orig + t * r = b.
		//  orign + t * r = 
		double a01 = b[0] - a[0];
		double a11 = b[1] - a[1];
		double det = r[0] * a11 - a01 * r[1];
		if (std::abs(det) > repsilon_) {
			//Noncolinear
			double s0 = b[0] - orig[0];
			double s1 = b[1] - orig[1];
			double t0 = (a11*s0 - a01 * s1) / det;
			double t1 = (r[0]*s1 - r[1] * s0) / det;
			return { t0,t1 };
		}
		else {
			//Colinear
			double t1;
			if (a01 > repsilon_) {
				t1 = (b[0] - orig[0]) / a01;
			}
			else {
				t1 = (b[1] - orig[1]) / a11;
			}
			if ((-repsilon_ < t1) & (t1 < 1 + repsilon_)) {
				return { 0.0, t1 };
			}
			else {
				if (r[0] > repsilon_) {
					return { (b[0] - orig[0]) / r[0], 0.0 };
				}
				else {

					return { (b[1] - orig[1]) / r[1], 0.0 };
				}
			}
		}
	}

	bool is_ray_segment_intersect(std::array<double, 2>orig, std::array<double, 2>r, std::array<double, 2>a, std::array<double, 2>b) {
		std::array<double, 2> t = ray_segment_intersect(orig, r, a, b);
		return (t[0] > -repsilon_) & (-repsilon_ < t[1]) & (t[1] < 1 + repsilon_);
	}

	bool is_point_in(std::array<double,2> point) {
		//Check by parity of a random ray trace
		double theta = M_PI * drand();
		std::array<double, 2> direct{cos(theta), sin(theta)};
		int isect_count = 0;
		for (int foo = 0; foo < boundary_->GetNumberOfCells(); foo++) {
			auto pa = boundary_->GetPoint(boundary_->GetCell(foo)->GetPointId(0));
			std::array<double, 2> paa {pa[0],pa[1]};
			auto pb = boundary_->GetPoint(boundary_->GetCell(foo)->GetPointId(1));
			std::array<double, 2> pba{ pb[0],pb[1] };
			std::array<double, 2> isect_params = ray_segment_intersect(point, direct, paa, pba);
			if ((isect_params[1] == 0.0)&(isect_params[0]>0.0)) {
				//This is slightly dangerous! Technically could be an infinite loop...
				//If the segment fully lies inside the ray, we reroll the ray direction and restart.
				std::cout << "Segment is subset of ray!";
				theta = M_PI * drand();
				direct={ cos(theta), sin(theta) };
				isect_count = 0;
				foo = -1;
				continue;
			}
			//Normally the ray hits the segment interior.
			if ((-repsilon_ < isect_params[1]) & (isect_params[1] < 1 + repsilon_)) {
				if (std::abs(isect_params[0]) <= repsilon_) {
					return true;
				}
				else if (isect_params[0] >0 ) {
					isect_count++;
				}
			}
		}
		return (isect_count % 2)==1;
	}

	double drand() {
		return 2.0*(((double) std::rand()) / RAND_MAX) - 1.0 ;
	}

	void test_rays() {
		std::ofstream ray_test_file;
		ray_test_file.open("raystest.txt");
		for (int foo = 0; foo < 42; foo++) {
			std::vector<double> somerands{ drand(), drand(), drand(), drand() , drand(), drand(), drand(), drand() };
			std::array<double, 2> orig{ somerands[0], somerands[1] };
			std::array<double, 2> r{ somerands[2], somerands[3] };
			std::array<double, 2> a{ somerands[4], somerands[5] };
			std::array<double, 2> b{ somerands[6], somerands[7] };
			auto t = ray_segment_intersect(orig, r, a, b);
			bool is_isect = is_ray_segment_intersect(orig, r, a, b);
			for (int bar = 0; bar < 8; bar++) ray_test_file << somerands[bar] << ", ";
			ray_test_file << is_isect << "\n";
		}
	}

	vtkSmartPointer<vtkDoubleArray> get_planar_normals() {
		vtkSmartPointer<vtkDoubleArray>  raw_normals = get_raw_normals();
		planar_normals_->SetNumberOfComponents(2);
		planar_normals_->SetNumberOfTuples(boundary_->GetNumberOfPoints());
		planar_normals_->SetName("Planar Normals");
		// Normalize raw normals if over a
		for (int foo = 0; foo < boundary_->GetNumberOfPoints(); foo++) {
			auto normal_foo = raw_normals->GetTuple(foo);
			planar_normals_->SetTuple(foo, normal_foo);
		}
		return planar_normals_;
	}

	vtkSmartPointer<vtkPolyData> get_boundary() {
		return boundary_;
	}
};

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
		double _corner[3], maxdir[3], middir[3], mindir[3], _sizes[3];
		auto obber = vtkSmartPointer<vtkOBBTree>::New();
		obber->ComputeOBB(above, _corner, maxdir, middir, mindir, _sizes);
		above_box_size_[0] = vtkMath::Norm(maxdir);
		above_box_size_[1] = vtkMath::Norm(middir);
		above_box_size_[2] = vtkMath::Norm(mindir);


		auto below = clipper->GetClippedOutput();
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
	vtkSmartPointer<vtkPolyData> plane_disks_ = vtkSmartPointer<vtkPolyData>::New();

public:
	vtkSmartPointer<vtkPolyData> get_plane_disks() {
		return plane_disks_;
	}

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
		//auto ppoints = vtkSmartPointer<vtkPoints2D>::New();
		//Project all points of the cross section polys into the plane with the {planarx_, planary_} orthonormal basis
		vtkMath::Subtract(disks->GetPoint(0), cut_origin_, planarx_);
		vtkMath::Normalize(planarx_);
		vtkMath::Cross(planenormal, planarx_, planary_);
		vtkMath::Normalize(planary_);
		//The ppoints are the in-plane points: planar points.
		auto pprojection = vtkSmartPointer<vtkTransform>::New();
		pprojection->PostMultiply();
		pprojection->Translate(-cut_origin_[0], -cut_origin_[1], -cut_origin_[2]);
		//const double ortho[16]{ 1.,0.,0.,0., 0.,1.,0.,0., 0.,0.,1.,0., 0.,0.,0.,1.};
		//const double ortho[16] {planarx_[0], planary_[0], planenormal[0], 0.0, planarx_[1], planary_[1], planenormal[1], 0.0, planarx_[2], planary_[2], planenormal[2], 0.0, 0.0, 0.0, 0.0, 1.0};
		const double ortho[16]{ planarx_[0],planarx_[1],planarx_[2],0.0, planary_[0],planary_[1],planary_[2],0.0, planenormal[0],planenormal[1],planenormal[2],0.0,  0.0,0.0,0.0,1.0};
		pprojection->Concatenate(ortho);
		auto transform_filter = vtkSmartPointer<vtkTransformPolyDataFilter>::New();
		transform_filter->SetTransform(pprojection);
		transform_filter->SetInputData(disks);
		transform_filter->Update();
		plane_disks_ = transform_filter->GetOutput();
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

void slicey_example() {
	boost::filesystem::path infilepath{ "C:\\Users\\sscott\\Pictures\\cube.stl" };
	auto reader = vtkSmartPointer<vtkSTLReader>::New();
	reader->SetFileName(infilepath.string().c_str());
	reader->Update();
	auto normaler = vtkSmartPointer<vtkPolyDataNormals>::New();
	normaler->SetInputConnection(reader->GetOutputPort());
	normaler->ComputeCellNormalsOn();
	normaler->ComputePointNormalsOn();
	normaler->Update();
	auto surface = normaler->GetOutput();



	auto sharpplane = vtkSmartPointer<vtkPlane>::New();
	//double anorigin[3]{ 100.0, 100.0, 10.0 };
	double anorigin[3]{ 0.100, 0.100, 0.100 };
	sharpplane->SetOrigin(anorigin);
	double adirection[3]{ 0.1235, 0.1235, 0.1235 };
	sharpplane->SetNormal(adirection);

	auto seamly = Seam(surface, sharpplane);
	auto curves = seamly.get_curves();

	auto writer = vtkSmartPointer<vtkXMLPolyDataWriter>::New();
	writer->SetInputData(curves);
	writer->SetFileName("aquickcurve.vtp");
	writer->Write();

	auto writer2 = vtkSmartPointer<vtkXMLPolyDataWriter>::New();
	writer2->SetInputData(seamly.get_disks());
	writer2->SetFileName("aquickdisk.vtp");
	writer2->Write();
}

vtkSmartPointer<vtkPolyData> generate_ngon(const int nn) {
	auto ngon_polydata = vtkSmartPointer<vtkPolyData>::New();
	auto ngon_pts = vtkSmartPointer<vtkPoints>::New();
	std::vector<std::vector<int>> cells;
	for (int foo = 0; foo < nn; foo++) {
		double ptcoor[3]{ cos(2 * M_PI / nn * foo), sin(2 * M_PI / nn * foo), 0.0 };
		ngon_pts->InsertPoint(foo, ptcoor);
		//std::vector<int> cell
		if ((foo > 0)&(foo < (nn - 1))) cells.push_back({ foo, foo + 1, 0 });
	}
	auto ngon_cells = vtkSmartPointer<vtkCellArray>::New();
	// Copy the triangles to the polydata
	for (auto cell : cells) {
		auto tri = vtkSmartPointer<vtkTriangle>::New();
		for (int foo = 0; foo < 3; foo++) tri->GetPointIds()->SetId(foo, cell[foo]);
		ngon_cells->InsertNextCell(tri);
	}
	ngon_polydata->SetPoints(ngon_pts);
	ngon_polydata->SetPolys(ngon_cells);
	return ngon_polydata;
}

vtkSmartPointer<vtkPolyData> generate_ngon_centered(const int nn) {
	auto ngon_polydata = vtkSmartPointer<vtkPolyData>::New();
	auto ngon_pts = vtkSmartPointer<vtkPoints>::New();
	std::vector<std::vector<int>> cells;
	double origin[3]{ 0.0, 0.0, 0.0 };
	ngon_pts->InsertPoint(nn, origin);
	for (int foo = 0; foo < nn; foo++) {
		double ptcoor[3]{ cos(2 * M_PI / nn * foo), sin(2 * M_PI / nn * foo), 0.0 };
		ngon_pts->InsertPoint(foo, ptcoor);
		//std::vector<int> cell
		cells.push_back({ foo, (foo + 1)%nn,  nn});
	}
	auto ngon_cells = vtkSmartPointer<vtkCellArray>::New();
	// Copy the triangles to the polydata
	for (auto cell : cells) {
		auto tri = vtkSmartPointer<vtkTriangle>::New();
		for (int foo = 0; foo < 3; foo++) tri->GetPointIds()->SetId(foo, cell[foo]);
		ngon_cells->InsertNextCell(tri);
	}
	ngon_polydata->SetPoints(ngon_pts);
	ngon_polydata->SetPolys(ngon_cells);
	return ngon_polydata;
}

vtkSmartPointer<vtkPolyData> generate_squarish(const int nn) {
	auto ngon_polydata = vtkSmartPointer<vtkPolyData>::New();
	auto ngon_pts = vtkSmartPointer<vtkPoints>::New();
	std::vector<std::vector<int>> cells;
	double origin[3]{ 0.0,0.0,0.0 };
	ngon_pts->InsertPoint(nn, origin);
	for (int foo = 0; foo < nn; foo++) {
		double ptcoor[3]{ cos(2 * M_PI / nn * foo), sin(2 * M_PI / nn * foo), 0.0 };
		vtkMath::MultiplyScalar2D(ptcoor, 1/(std::abs(ptcoor[0]) + std::abs(ptcoor[1])) );
		ngon_pts->InsertPoint(foo, ptcoor);
		//std::vector<int> cell
		cells.push_back({ foo, (foo + 1) % nn,  nn });
	}
	auto ngon_cells = vtkSmartPointer<vtkCellArray>::New();
	// Copy the triangles to the polydata
	for (auto cell : cells) {
		auto tri = vtkSmartPointer<vtkTriangle>::New();
		for (int foo = 0; foo < 3; foo++) tri->GetPointIds()->SetId(foo, cell[foo]);
		ngon_cells->InsertNextCell(tri);
	}
	ngon_polydata->SetPoints(ngon_pts);
	ngon_polydata->SetPolys(ngon_cells);
	return ngon_polydata;
}

vtkSmartPointer<vtkPolyData> generate_star(const int num_outer_points = 5, const double in_radius_multiplier = -1.0) {
	std::vector<ppoint> star_ppts;
	star_ppts.resize(2 * num_outer_points);
	for (int foo = 0; foo < num_outer_points; foo++) {
		star_ppts[foo] = { -sin(2 * M_PI / num_outer_points * foo), cos(2 * M_PI / num_outer_points * foo) };
	}
	auto star_polydata = vtkSmartPointer<vtkPolyData>::New();
	PlanarNormalFilter pnf = PlanarNormalFilter(star_polydata);
	for (int foo = 0; foo < num_outer_points; foo++) {
		ppoint orig = star_ppts[foo];
		ppoint dir = star_ppts[(foo+2) % num_outer_points] - star_ppts[foo];
		ppoint a = star_ppts[(foo + num_outer_points- 1) % num_outer_points];
		ppoint b = star_ppts[(foo + 1) % num_outer_points];
		ppoint st = pnf.ray_segment_intersect(orig, dir, a, b);
		star_ppts[foo + num_outer_points] = orig + st[0] * dir;
	}
	if (in_radius_multiplier > 0.0) {
		for (int foo = 0; foo < num_outer_points; foo++) {
			star_ppts[foo + num_outer_points] = in_radius_multiplier * star_ppts[foo + num_outer_points];
		}
	}
	auto star_polydata_pts = vtkSmartPointer<vtkPoints>::New();
	for (int foo = 0; foo < 2*num_outer_points; foo++) {
		double pfoo[3]{ star_ppts[foo][0], star_ppts[foo][1], 0.0 };
		star_polydata_pts->InsertNextPoint(pfoo);
	}
	auto star_cells = vtkSmartPointer<vtkCellArray>::New();
	// Copy the edges to the polydata
	for (int foo = 0; foo < num_outer_points; foo++) {
		auto edge = vtkSmartPointer<vtkLine>::New();
		edge->GetPointIds()->SetId(0, foo);
		edge->GetPointIds()->SetId(1, foo + num_outer_points);
		star_cells->InsertNextCell(edge);
		auto edge2 = vtkSmartPointer<vtkLine>::New();
		edge2->GetPointIds()->SetId(0, foo);
		edge2->GetPointIds()->SetId(1, (foo+ num_outer_points -1) % num_outer_points + num_outer_points);
		star_cells->InsertNextCell(edge2);
	}
	star_polydata->SetPoints(star_polydata_pts);
	star_polydata->SetLines(star_cells);
	return star_polydata;
}

vtkSmartPointer<vtkPolyData> subdivide_edges(const vtkSmartPointer<vtkPolyData> x, const int splits = 2) {
	auto y = vtkSmartPointer<vtkPolyData>::New();
	auto y_pts = vtkSmartPointer<vtkPoints>::New();
	int counter_pnts = 0;
	for (int foo = 0; foo < x->GetNumberOfPoints(); foo++) {
		y_pts->InsertNextPoint(x->GetPoint(foo));
		counter_pnts++;
	}
	auto y_lines = vtkSmartPointer<vtkCellArray>::New();
	auto nlines = x->GetNumberOfLines();
	for (int foo = 0; foo < x->GetNumberOfLines(); foo++) {
		double pta[3], ptb[3];
		auto ptids = vtkSmartPointer<vtkIdList>::New();
		std::cout << "cell " << foo << "\n";
		x->GetLines()->GetCellAtId(foo, ptids);
		int id_a = ptids->GetId(0);
		x->GetPoints()->GetPoint(id_a, pta);
		int id_b = ptids->GetId(1);
		x->GetPoint(id_b, ptb);
		//std::vector<ppoint> newpoints;
		//newpoints.resize(splits + 1);
		ppoint a{ pta[0],pta[1] };
		ppoint b{ ptb[0], ptb[1] };
		//newpoints[0] = a;
		//newpoints[splits] = b;
		std::vector<int> pt_order;
		pt_order.push_back(id_a);
		for (int bar = 1; bar < splits; bar++) {
			double t = ((double)bar) / splits;
			ppoint newp = t * b + (1 - t) * a;
			pt_order.push_back(counter_pnts);
			y_pts->InsertNextPoint(newp[0],newp[1], 0.0);
			counter_pnts++;
			pt_order.push_back(counter_pnts);
		}
		pt_order.push_back(id_b);
		for (int bar = 0; bar < splits; bar++) {
			auto edge = vtkSmartPointer<vtkLine>::New();
			edge->GetPointIds()->SetId(0, pt_order[bar]);
			edge->GetPointIds()->SetId(1, pt_order[bar+1]);
			y_lines->InsertNextCell(edge);
		}
	}
	y->SetPoints(y_pts);
	y->SetLines(y_lines);
	return y;
}

void write_it(const vtkSmartPointer<vtkPolyData> x, const std::string filename) {
	//Write
	auto writer3 = vtkSmartPointer<vtkXMLPolyDataWriter>::New();
	writer3->SetInputData(x);
	writer3->SetFileName(filename.c_str());
	writer3->Write();
}

void plane_normals_example(const int ngonnum) {	
	auto ngon = generate_ngon(ngonnum);
	PlanarNormalFilter normalizer = PlanarNormalFilter(ngon);
	vtkSmartPointer<vtkDoubleArray> normal_section = normalizer.get_planar_normals();

	vtkSmartPointer<vtkPolyData> bounding_curve = normalizer.get_boundary();
	bounding_curve->GetPointData()->AddArray(normal_section);
	write_it(bounding_curve, "bounding_curve.vtp");

	ngon->GetPointData()->AddArray(normal_section);
	write_it(ngon, "ngon.vtp");
}

void example_ngon_normals2(const int ngonnum) {
	auto ngon = generate_ngon_centered(ngonnum);
	PlanarNormalFilter normalizer = PlanarNormalFilter(ngon);
	vtkSmartPointer<vtkDoubleArray> normal_section = normalizer.get_planar_normals();

	vtkSmartPointer<vtkPolyData> bounding_curve = normalizer.get_boundary();
	bounding_curve->GetPointData()->AddArray(normal_section);
	write_it(bounding_curve, "bounding_curve.vtp");

	ngon->GetPointData()->AddArray(normal_section);
	write_it(ngon, "ngon.vtp");
}

void example_square_normals(const int ngonnum) {
	auto ngon = generate_squarish(ngonnum);
	PlanarNormalFilter normalizer = PlanarNormalFilter(ngon);
	vtkSmartPointer<vtkDoubleArray> normal_section = normalizer.get_planar_normals();

	vtkSmartPointer<vtkPolyData> bounding_curve = normalizer.get_boundary();
	bounding_curve->GetPointData()->AddArray(normal_section);
	write_it(bounding_curve, "bounding_curve.vtp");
}

void example_ray_segger() {
	auto ngon = generate_ngon(3);
	PlanarNormalFilter normalizer = PlanarNormalFilter(ngon);
	normalizer.test_rays();
}

void example_in_point_check0() {
	boost::filesystem::path infilepath{ "C:\\Users\\sscott\\Pictures\\trex_connected.stl" };
	auto reader = vtkSmartPointer<vtkSTLReader>::New();
	reader->SetFileName(infilepath.string().c_str());
	reader->Update();
	auto surface = reader->GetOutput();
	auto sharpplane = vtkSmartPointer<vtkPlane>::New();
	double anorigin[3]{ 100.0, 100.0, 10.0 };
	//double anorigin[3]{ 0.100, 0.100, 0.100 };
	sharpplane->SetOrigin(anorigin);
	double adirection[3]{ 0.1235, 0.1235, 0.1235 };
	sharpplane->SetNormal(adirection);

	auto seamly = Seam(surface, sharpplane);
	auto cross_section = CrossSection(sharpplane, seamly.get_disks());
	auto planar_filter = PlanarNormalFilter(cross_section.get_plane_disks());
	auto curve = planar_filter.get_boundary();

	std::ofstream p_file;
	p_file.open("point_test_fielda.txt");

	for (int foo = 0; foo < curve->GetNumberOfPoints(); foo++) {
		auto pt_i = curve->GetPoint(foo);
		std::array<double, 2> pt {pt_i[0], pt_i[1]};
		bool is_in = planar_filter.is_point_in(pt);
		p_file << pt[0] << "," << pt[1] << "," << 0.0 << "," << is_in << "\n";
	}
	const int num_added_points = 3000;
	const int num_points_started = curve->GetNumberOfPoints();
	for (int foo = num_points_started; foo < num_points_started + num_added_points; foo++) {
		//double pt_i[3]{ 0.208747, 3.57491, 0.0 };
		double pt_i[3] {8*planar_filter.drand(), 11*planar_filter.drand(), 0 };
		std::array<double, 2> pt{ pt_i[0], pt_i[1] };
		bool is_in = planar_filter.is_point_in(pt);
		p_file << pt[0] << "," << pt[1] << "," << 0.0 << "," << is_in << "\n";
	}
	write_it(curve, "example_in_point_check0.vtp");
}

void example_in_point_check_on_edge() {
	auto ngon = generate_ngon(3);
	auto planar_filter = PlanarNormalFilter(ngon);
	auto curve = planar_filter.get_boundary();

	std::ofstream p_file;
	p_file.open("point_test_edge.txt");
 
	const int num_pt_per_edge = 360;

	for (int foo = 0; foo < curve->GetNumberOfCells(); foo++) {
	}

	for (int foo = 0; foo < curve->GetNumberOfCells(); foo++) {
		auto pts_in_i = curve->GetCell(foo)->GetPointIds();
		double pa[3], pb[3];
		curve->GetPoint(pts_in_i->GetId(0), pa);
		curve->GetPoint(pts_in_i->GetId(1), pb);
		for (int foo = 0; foo < num_pt_per_edge; foo++) {
			double tt = std::abs(planar_filter.drand());
			double pt_i[3]{ 0.0,0.0,0.0};
			pt_i[0] = pa[0] * tt + pb[0] * (1 - tt);
			pt_i[1] = pa[1] * tt + pb[1] * (1 - tt);
			std::array<double, 2> pt{ pt_i[0], pt_i[1] };
			bool is_in = planar_filter.is_point_in(pt);
			p_file << pt[0] << "," << pt[1] << "," << 0.0 << "," << is_in << "\n";
		}
	}

	write_it(curve, "example_in_point_check_on_edge.vtp");
}


void example_star_maker(const int starpointed=5, const double in_radius_multipier=-1.0, const int splits = 2) {
	auto star = generate_star(starpointed, in_radius_multipier);
	//auto starlines = star->GetLines();
	//for (int foo = 0; foo < star->GetNumberOfLines(); foo++){
	//	std::cout << foo << "\n";
	//	auto id_list = vtkSmartPointer<vtkIdList>::New();
	//	starlines->GetCellAtId(foo, id_list);

	//}
	auto star_sub = subdivide_edges(star, splits);
	write_it(star_sub, "star.vtp");
}

int main() {
	example_star_maker(8,0.5,3);
};