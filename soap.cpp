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
#include <vtkAppendPolyData.h>
#include <vtkCenterOfMass.h>
#include <vtkCellArray.h>
#include <vtkCellLocator.h>
#include <vtkCleanPolyData.h>
#include <vtkClipPolyData.h>
#include <vtkContourFilter.h>
#include <vtkContourTriangulator.h>
#include <vtkCutter.h>
#include <vtkDelaunay2D.h>
#include <vtkDistancePolyDataFilter.h>
#include <vtkDoubleArray.h>
#include <vtkFeatureEdges.h>
#include <vtkImplicitPolyDataDistance.h>
#include <vtkLine.h>
#include <vtkLinearExtrusionFilter.h>
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

ostream & operator << (ostream &out, const ppoint p) {
	out << p[0] << " " << p[1];
	return out;
}

void write_it(const vtkSmartPointer<vtkPolyData> x, const std::string filename) {
	//Write
	auto writer3 = vtkSmartPointer<vtkXMLPolyDataWriter>::New();
	writer3->SetInputData(x);
	writer3->SetFileName(filename.c_str());
	writer3->Write();
}


class PlanarNormalFilter {
private:
	vtkSmartPointer<vtkPolyData> polys_ = vtkSmartPointer<vtkPolyData>::New();
	vtkSmartPointer<vtkPolyData> boundary_ = vtkSmartPointer<vtkPolyData>::New();
	vtkSmartPointer<vtkDoubleArray> planar_normals_ = vtkSmartPointer<vtkDoubleArray>::New();
	double min_boundary_edge_length_{ -1.0 };
	const double repsilon_ = std::sqrt(std::numeric_limits<double>::epsilon());

public:
	PlanarNormalFilter() {
	};

	PlanarNormalFilter(vtkSmartPointer<vtkPolyData> polys) {
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
				//Colinear
			normal_i[0] = pa[1] - pc[1];
			normal_i[1] = pc[0] - pa[0];
			normal_i[2] = 0.0;
			vtkMath::Normalize2D(normal_i);
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
		double push_off_dist = std::max(min_boundary_edge_length_/2, 100*repsilon_);
		for (int foo = 0; foo < boundary_->GetNumberOfPoints(); foo++) {
			double normal_foo[2], atpnt[3];
			raw_normals->GetTuple(foo, normal_foo);
			boundary_->GetPoint(foo, atpnt);
			ppoint lil_off = { atpnt[0] + push_off_dist * normal_foo[0], atpnt[1] + push_off_dist * normal_foo[1] };
			bool is_in = is_point_in(lil_off);
			if (is_in) {
				double neg_normal[3]{ -normal_foo[0], -normal_foo[1], 0.0};
				planar_normals_->SetTuple(foo, neg_normal);
			}
			else {
				planar_normals_->SetTuple(foo, normal_foo);
			}
		}
		return planar_normals_;
	}

	void set_boundary(vtkSmartPointer<vtkPolyData> boundary) {
		boundary_ = boundary;
	}

	vtkSmartPointer<vtkPolyData> get_boundary() {
		return boundary_;
	}

	vtkSmartPointer<vtkPolyData> offsetter_field(double min_dist, double max_dist, int sampling = 23) {
		vtkSmartPointer<vtkDoubleArray> normals = get_planar_normals();
		auto off_pnts = vtkSmartPointer<vtkPoints>::New();
		for (int foo = 0; foo < boundary_->GetNumberOfPoints(); foo++) {
			double pnt_foo[3];
			boundary_->GetPoint(foo, pnt_foo);
			off_pnts->InsertNextPoint(pnt_foo);
			ppoint at{ pnt_foo[0], pnt_foo[1] };
			//
			double nrm_foo[2];
			normals->GetTuple(foo, nrm_foo);
			ppoint nrm { nrm_foo[0], nrm_foo[1] };
			for (int bar = 0; bar < sampling; bar++) {
				double ss = bar / (sampling - 1.0);
				double tt = max_dist * ss + min_dist * (1 - ss);
				ppoint nat = at + tt * nrm;
				off_pnts->InsertNextPoint(nat[0], nat[1], 0.0);
			}
		}
		auto delaunay = vtkSmartPointer<vtkDelaunay2D>::New();
		auto just_points = vtkSmartPointer<vtkPolyData>::New();
		just_points->SetPoints(off_pnts);
		delaunay->SetInputData(just_points);
		delaunay->Update();
		auto offgrid = delaunay->GetOutput();

		//boundary_->SetCells(cellarray)
		//return delaunay->GetOutput();
		auto dister = vtkSmartPointer<vtkImplicitPolyDataDistance>::New();
		dister->SetInput(boundary_);
		
		auto dist = vtkSmartPointer<vtkDoubleArray>::New();
		dist->SetNumberOfComponents(1);
		dist->SetName("offset");
		auto cell_locater = vtkSmartPointer<vtkCellLocator>::New();
		cell_locater->SetDataSet(boundary_);
		cell_locater->BuildLocator();
		for (int foo = 0; foo < offgrid->GetNumberOfPoints(); foo++) {
			double pnt_foo[3];
			offgrid->GetPoint(foo, pnt_foo);
			double closest_point[3], dist2;
			vtkIdType cellid;
			int subid;
			cell_locater->FindClosestPoint(pnt_foo, closest_point, cellid, subid, dist2);
			bool is_in = is_point_in({ pnt_foo[0], pnt_foo[1] });
			if (is_in) {
				dist->InsertNextValue(-std::sqrt(dist2));
			}
			else {
				dist->InsertNextValue(std::sqrt(dist2));
			}
		}
		offgrid->GetPointData()->AddArray(dist);
		return offgrid;
	}
};

struct SeamParameters {
	double groove_outer = 1.0;
	double groove_inner = - 1.0;
	double gap_radial = 0.5;
	double gap_depth = 0.5;
	double tongue_depth = 5.0;
	double trim_depth = 1.0;
};

class Expandilizer {
private:
	vtkSmartPointer<vtkPolyData> offgrid_ = vtkSmartPointer<vtkPolyData>::New();
	vtkSmartPointer<vtkContourFilter> isoer_ = vtkSmartPointer<vtkContourFilter>::New();

public:
	Expandilizer(vtkSmartPointer<vtkPolyData> offgrid) {
		offgrid_ = offgrid;
		write_it(offgrid_, "offgrid_was.vtp");
		offgrid_->GetPointData()->SetActiveScalars("offset");
		isoer_->SetInputData(offgrid_);
		isoer_->ComputeScalarsOff();
	}

	vtkSmartPointer<vtkPolyData> flat(double off, double z) {
		isoer_->SetNumberOfContours(1);
		isoer_->SetValue(0, off);
		isoer_->Update();
		auto filler = vtkSmartPointer<vtkContourTriangulator>::New();
		filler->SetInputConnection(isoer_->GetOutputPort());
		filler->Update();
		auto pprojection = vtkSmartPointer<vtkTransform>::New();
		pprojection->Translate(0.0, 0.0, z);
		auto transform_filter = vtkSmartPointer<vtkTransformPolyDataFilter>::New();
		transform_filter->SetTransform(pprojection);
		transform_filter->SetInputData(filler->GetOutput());
		transform_filter->Update();
		return transform_filter->GetOutput();
	}

	vtkSmartPointer<vtkPolyData> flat(double off1, double off2, double z) {
		isoer_->SetNumberOfContours(2);
		isoer_->SetValue(0, off1);
		isoer_->SetValue(1, off2);
		isoer_->Update();
		auto filler = vtkSmartPointer<vtkContourTriangulator>::New();
		filler->SetInputConnection(isoer_->GetOutputPort());
		filler->Update();
		auto pprojection = vtkSmartPointer<vtkTransform>::New();
		pprojection->Translate(0.0, 0.0, z);
		auto transform_filter = vtkSmartPointer<vtkTransformPolyDataFilter>::New();
		transform_filter->SetTransform(pprojection);
		transform_filter->SetInputData(filler->GetOutput());
		transform_filter->Update();
		return transform_filter->GetOutput();
	}

	vtkSmartPointer<vtkPolyData> band(double off, double z1, double z2) {
		isoer_->SetNumberOfContours(1);
		isoer_->SetValue(0, off);
		isoer_->Update();
		auto extruder = vtkSmartPointer<vtkLinearExtrusionFilter>::New();
		extruder->SetInputConnection(isoer_->GetOutputPort());
		extruder->SetVector(0.0, 0.0, 1.0);
		extruder->SetScaleFactor(z2 - z1);
		extruder->Update();
		auto pprojection = vtkSmartPointer<vtkTransform>::New();
		pprojection->Translate(0.0, 0.0, z1);
		auto transform_filter = vtkSmartPointer<vtkTransformPolyDataFilter>::New();
		transform_filter->SetTransform(pprojection);
		transform_filter->SetInputData(extruder->GetOutput());
		transform_filter->Update();
		return transform_filter->GetOutput();
	}

	vtkSmartPointer<vtkPolyData> tongue_and_groove(SeamParameters params) {
		auto appender = vtkSmartPointer<vtkAppendPolyData>::New();
		//The isocontour extraction is super redundant here fix in future!!!!
		double z_tongue_tip = (params.gap_depth - params.tongue_depth) / 2.0;
		double z_teeth_back = z_tongue_tip - params.gap_depth;
		double z_teeth_front = z_teeth_back - params.trim_depth;
		auto flat1 = flat(params.groove_outer, z_teeth_front);
		appender->AddInputData(flat1);
		auto band2 = band(params.groove_outer, z_teeth_front, -z_tongue_tip);
		appender->AddInputData(band2);
		auto flat3 = flat(params.groove_inner, params.groove_outer, -z_tongue_tip);
		appender->AddInputData(flat3);
		auto band4 = band(params.groove_inner, z_teeth_back, -z_tongue_tip);
		appender->AddInputData(band4);
		auto flat5 = flat(params.groove_inner, z_teeth_back);
		appender->AddInputData(flat5);
		double innest = params.groove_inner - params.gap_radial;
		auto flat6 = flat(innest, z_tongue_tip);
		appender->AddInputData(flat6);
		auto band7 = band(innest, z_tongue_tip, -z_teeth_back);
		appender->AddInputData(band7);
		auto flat8 = flat(innest, params.groove_outer, -z_teeth_back);
		appender->AddInputData(flat8);
		auto band9 = band(params.groove_outer, -z_teeth_back, -z_teeth_front);
		appender->AddInputData(band9);
		auto flat10 = flat(params.groove_outer, -z_teeth_front);
		appender->AddInputData(flat10);
		appender->Update();
		auto cleaner = vtkSmartPointer<vtkCleanPolyData>::New();
		cleaner->PointMergingOn();
		cleaner->SetInputData(appender->GetOutput());
		cleaner->Update();
		return cleaner->GetOutput();
		//Put all that together into one mesh and merge all the duplicated points...
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
	vtkSmartPointer<vtkTransform> transform_ = vtkSmartPointer<vtkTransform>::New();
	std::array<double, 3> tongue_direction_{ 0.0, 0.0, 0.0 };

public:
	vtkSmartPointer<vtkPolyData> get_plane_disks() {
		return plane_disks_;
	}
	
	vtkSmartPointer<vtkTransform> get_transform() {
		return transform_;
	}

	std::array<double, 3> get_cut_origin() {
		return {cut_origin_[0], cut_origin_[1], cut_origin_[2]};
	}


	vtkSmartPointer<vtkTransform> get_sheer_plane_transform(std::array<double,3> cut_direction, std::array<double,3> tongue_direction) {
		// If the slide direction is not the cut plane direction, then you want to sheer the cross section into
		// a plane whose normal is slide direction, design the joint, the resheer it back in the cut-plane.
		double adirection[3]{cut_direction[0], cut_direction[1], cut_direction[2]};
		double sdirection[3]{ tongue_direction[0], tongue_direction[1], tongue_direction[2] };
		double dircross[3];
		vtkMath::Normalize(sdirection);
		vtkMath::Normalize(adirection);
		vtkMath::Cross(sdirection, adirection, dircross);
		double sin_dir_angle = vtkMath::Norm(dircross);
		double cos_dir_angle = vtkMath::Dot(adirection, sdirection);
		vtkMath::Normalize(dircross);
		double dir_basis[3][3]{ {sdirection[0], adirection[0], dircross[0]}, {sdirection[1], adirection[1], dircross[1]}, {sdirection[2], adirection[2], dircross[2]} };
		double dir_basis_inv[3][3];
		vtkMath::Invert3x3(dir_basis, dir_basis_inv);
		double tar_basis[3][3]{ {0.0, sin_dir_angle, 0.0},{0.0,0.0,1.0},{1.0,cos_dir_angle,0.0} };
		double qmat[3][3], qmat_inv[3][3], sheer_mat[3][3], qmat_tmp[3][3];
		double sheer_mat_lower[3][3]{ {1.0, 0.0, 0.0},{0.0, 1.0, 0.0}, {sin_dir_angle/cos_dir_angle, 0.0, 1.0} };
		vtkMath::Multiply3x3(tar_basis, dir_basis_inv, qmat_tmp);
		vtkMath::Orthogonalize3x3(qmat_tmp, qmat);
		vtkMath::Transpose3x3(qmat, qmat_inv);
		double tempmat[3][3];
		vtkMath::Multiply3x3(qmat_inv, sheer_mat_lower, tempmat);
		vtkMath::Multiply3x3(tempmat, qmat, sheer_mat);
		auto smosher = vtkSmartPointer<vtkTransform>::New();
		smosher->PostMultiply();
		smosher->Translate(-cut_origin_[0], -cut_origin_[1], -cut_origin_[2]);
		double sheer_mat_elements[16]{ sheer_mat[0][0], sheer_mat[0][1], sheer_mat[0][2], 0.0, sheer_mat[1][0], sheer_mat[1][1], sheer_mat[1][2], 0.0, sheer_mat[2][0], sheer_mat[2][1], sheer_mat[2][2], 0.0, 0.0, 0.0, 0.0,1.0 };
		smosher->Concatenate(sheer_mat_elements);
		smosher->Translate(cut_origin_[0], cut_origin_[1], cut_origin_[2]);
		return smosher;
	}

	vtkSmartPointer<vtkTransform> get_sheer_plane_transform() {
		double cut_plane_direction[3];
		plane_->GetNormal(cut_plane_direction);
		return get_sheer_plane_transform({ cut_plane_direction[0], cut_plane_direction[1] , cut_plane_direction[2] }, { tongue_direction_[0], tongue_direction_[1] , tongue_direction_[2] });
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
		transform_->PostMultiply();
		transform_->Translate(-cut_origin_[0], -cut_origin_[1], -cut_origin_[2]);
		//const double ortho[16]{ 1.,0.,0.,0., 0.,1.,0.,0., 0.,0.,1.,0., 0.,0.,0.,1.};
		//const double ortho[16] {planarx_[0], planary_[0], planenormal[0], 0.0, planarx_[1], planary_[1], planenormal[1], 0.0, planarx_[2], planary_[2], planenormal[2], 0.0, 0.0, 0.0, 0.0, 1.0};
		const double ortho[16]{ planarx_[0],planarx_[1],planarx_[2],0.0, planary_[0],planary_[1],planary_[2],0.0, planenormal[0],planenormal[1],planenormal[2],0.0,  0.0,0.0,0.0,1.0};
		transform_->Concatenate(ortho);
		auto transform_filter = vtkSmartPointer<vtkTransformPolyDataFilter>::New();
		transform_filter->SetTransform(transform_);
		transform_filter->SetInputData(disks);
		transform_filter->Update();
		plane_disks_ = transform_filter->GetOutput();
	}

	CrossSection(vtkSmartPointer<vtkPlane> plane, vtkSmartPointer<vtkPolyData> disks, std::array<double,3> tongue_direction) {
		plane_ = plane;
		tongue_direction_ = tongue_direction;
		vtkMath::Normalize(tongue_direction_.data());
		//Compute cut origin
		double origin_tmp0[3];
		vtkCenterOfMass::ComputeCenterOfMass(disks->GetPoints(), nullptr, origin_tmp0);
		plane_->ProjectPoint(origin_tmp0, cut_origin_);
		transform_ = get_sheer_plane_transform();		
		//The ppoints are the in-plane points: planar points.
		//auto ppoints = vtkSmartPointer<vtkPoints2D>::New();
		//Project all points of the cross section polys into the plane with the {planarx_, planary_} orthonormal basis
		double tmp_Tpoint[3];
		double tmp_point[3];
		disks->GetPoint(0, tmp_point);
		transform_->TransformPoint(tmp_point, tmp_Tpoint);
		vtkMath::Subtract(tmp_Tpoint, cut_origin_, planarx_);
		vtkMath::Normalize(planarx_);
		vtkMath::Cross(tongue_direction_.data(), planarx_, planary_);
		vtkMath::Normalize(planary_);
		//The ppoints are the in-plane points: planar points.
		transform_->PostMultiply();
		transform_->Translate(-cut_origin_[0], -cut_origin_[1], -cut_origin_[2]);
		const double ortho[16]{ planarx_[0],planarx_[1],planarx_[2],0.0, planary_[0],planary_[1],planary_[2],0.0, tongue_direction_[0], tongue_direction_[1], tongue_direction_[2],0.0,  0.0,0.0,0.0,1.0 };
		transform_->Concatenate(ortho);
		auto transform_filter = vtkSmartPointer<vtkTransformPolyDataFilter>::New();
		transform_filter->SetTransform(transform_);
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
		x->GetLines()->GetCellAtId(foo, ptids);
		int id_a = ptids->GetId(0);
		x->GetPoints()->GetPoint(id_a, pta);
		int id_b = ptids->GetId(1);
		x->GetPoint(id_b, ptb);
		ppoint a{ pta[0], pta[1] };
		ppoint b{ ptb[0], ptb[1] };
		std::vector<int> pt_order;
		pt_order.push_back(id_a);
		for (int bar = 1; bar < splits; bar++) {
			double t = ((double)bar) / splits;
			ppoint newp = t * b + (1 - t) * a;
			pt_order.push_back(counter_pnts);
			y_pts->InsertNextPoint(newp[0],newp[1], 0.0);
			counter_pnts++;
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

void example_plane_normals(const int ngonnum) {	
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

void example_in_point_check_star_nearly() {
	std::ofstream p_file;
	p_file.open("star_in_out_test_points.txt");
	const int starpointed = 7;
	const double in_radius_multipier = 0.7;
	const int splits = 33;
	auto star = generate_star(starpointed, in_radius_multipier);
	auto star_sub = subdivide_edges(star, splits);
	PlanarNormalFilter normalizer = PlanarNormalFilter();
	normalizer.set_boundary(star_sub);
	vtkSmartPointer<vtkDoubleArray> normal_section = normalizer.get_planar_normals();
	double repsilon = 100000*std::sqrt(std::numeric_limits<double>::epsilon());
	for (int foo = 0; foo < star_sub->GetNumberOfPoints(); foo++) {
		double normal_foo[2], atpnt[3];
		normal_section->GetTuple(foo, normal_foo);
		star_sub->GetPoint(foo, atpnt);
		for (int tt = -5; tt < 6; tt++) {
			ppoint lil_off = { atpnt[0] + tt * repsilon * normal_foo[0], atpnt[1] + tt * repsilon * normal_foo[1] };
			bool is_in = normalizer.is_point_in(lil_off);
			if (is_in) {
				p_file << lil_off[0] << ", " << lil_off[1] << ", " << 0.0 << ", " << 1 << "\n";
			}
			else {
				p_file << lil_off[0] << ", " << lil_off[1] << ", " << 0.0 << ", " << 0 << "\n";
			}
		}
	}
	p_file.close();
	star_sub->GetPointData()->AddArray(normal_section);
	write_it(star_sub, "star.vtp");
}

void example_star_maker(const int starpointed=5, const double in_radius_multipier=-1.0, const int splits = 2) {
	auto star = generate_star(starpointed, in_radius_multipier);
	auto star_sub = subdivide_edges(star, splits);
	PlanarNormalFilter normalizer = PlanarNormalFilter();
	normalizer.set_boundary(star_sub);
	vtkSmartPointer<vtkDoubleArray> normal_section = normalizer.get_planar_normals();
	star_sub->GetPointData()->AddArray(normal_section);
	write_it(star_sub, "star.vtp");
}

void example_offset_field() {
	auto star_simple = generate_star(7, 0.7);
	//Contrast with subdivide_edges(star_simple, 3);
	//Example shows the need to ensure edges are sufficiently subdivided.
	//Otherwise fake isolines can show up.
	auto star = subdivide_edges(star_simple, 33);
	write_it(star, "star.vtp");
	PlanarNormalFilter normalizer = PlanarNormalFilter();
	normalizer.set_boundary(star);
	vtkSmartPointer<vtkPolyData> off_field = normalizer.offsetter_field(-0.11, 0.11);
	write_it(off_field, "offstar.vtp");
}

void example_lattice_slicey() {
	boost::filesystem::path infilepath{ "C:\\Users\\sscott\\Pictures\\lattice1.stl" };
	auto reader = vtkSmartPointer<vtkSTLReader>::New();
	reader->SetFileName(infilepath.string().c_str());
	reader->Update();
	auto surface = reader->GetOutput();

	auto sharpplane = vtkSmartPointer<vtkPlane>::New();
	double anorigin[3]{684.4189787288672, 0.0, 0.0};
	//double anorigin[3]{ 100.0, 100.0, 10.0 };
	sharpplane->SetOrigin(anorigin);
	double adirection[3]{-0.5394782755183589, -0.7809840637246965, -0.3146856883491796};
	//double adirection[3]{ 0.1235, 0.1235, 0.1235 };
	sharpplane->SetNormal(adirection);

	auto seamly = Seam(surface, sharpplane);
	auto cross_sectly = CrossSection(sharpplane, seamly.get_disks());
	//write_it(cross_sectly.get_plane_disks(), "hipslice.vtp");
	PlanarNormalFilter planer = PlanarNormalFilter(cross_sectly.get_plane_disks());
	auto normals = planer.get_planar_normals();
	auto boundwithnormals = planer.get_boundary();
	boundwithnormals->GetPointData()->AddArray(normals);
	boundwithnormals->GetPointData()->AddArray(planer.get_raw_normals());
	write_it(boundwithnormals, "check_dez_normals.vtp");

	auto offler = planer.offsetter_field(-2.0, 5.0);
	write_it(offler, "hipofffield.vtp");
	auto expandilizer = Expandilizer(offler);
	SeamParameters params;
	auto widget = expandilizer.tongue_and_groove(params);
	write_it(widget, "widget.vtp");
}

void example_lattice_join_t() {
	boost::filesystem::path infilepath{ "C:\\Users\\sscott\\Pictures\\lattice1.stl" };
	auto reader = vtkSmartPointer<vtkSTLReader>::New();
	reader->SetFileName(infilepath.string().c_str());
	reader->Update();
	auto surface = reader->GetOutput();

	auto sharpplane = vtkSmartPointer<vtkPlane>::New();
	double anorigin[3]{ 684.4189787288672, 0.0, 0.0 };
	//double anorigin[3]{ 100.0, 100.0, 10.0 };
	sharpplane->SetOrigin(anorigin);
	double adirection[3]{ -0.5394782755183589, -0.7809840637246965, -0.3146856883491796 };
	//double adirection[3]{ 0.1235, 0.1235, 0.1235 };
	sharpplane->SetNormal(adirection);

	auto seamly = Seam(surface, sharpplane);
	auto cross_sectly = CrossSection(sharpplane, seamly.get_disks());
	//write_it(cross_sectly.get_plane_disks(), "hipslice.vtp");
	PlanarNormalFilter planer = PlanarNormalFilter(cross_sectly.get_plane_disks());
	auto normals = planer.get_planar_normals();
	auto boundwithnormals = planer.get_boundary();
	boundwithnormals->GetPointData()->AddArray(normals);
	boundwithnormals->GetPointData()->AddArray(planer.get_raw_normals());
	write_it(boundwithnormals, "check_dez_normals.vtp");

	auto offler = planer.offsetter_field(-2.0, 5.0);
	write_it(offler, "hipofffield.vtp");
	auto expandilizer = Expandilizer(offler);
	SeamParameters params;
	auto widget = expandilizer.tongue_and_groove(params);

	vtkSmartPointer<vtkTransform> transform = cross_sectly.get_transform();
	auto inv_transform = vtkSmartPointer<vtkTransform>::New();
	auto inv_transform_matrix = vtkSmartPointer<vtkMatrix4x4>::New();
	transform->GetInverse(inv_transform_matrix);
	inv_transform->SetMatrix(inv_transform_matrix);
	auto inver = vtkSmartPointer<vtkTransformPolyDataFilter>::New();
	inver->SetInputData(widget);
	inver->SetTransform(inv_transform);
	inver->Update();
	auto theseam = inver->GetOutput();
	write_it(theseam, "seam_seam_a_song.vtp");
}

void example_lattice_join_different_slide_direction() {
	boost::filesystem::path infilepath{ "C:\\Users\\sscott\\Pictures\\lattice1.stl" };
	auto reader = vtkSmartPointer<vtkSTLReader>::New();
	reader->SetFileName(infilepath.string().c_str());
	reader->Update();
	auto surface = reader->GetOutput();

	auto sharpplane = vtkSmartPointer<vtkPlane>::New();
	double anorigin[3]{ 684.4189787288672, 0.0, 0.0 };
	sharpplane->SetOrigin(anorigin);
	double adirection[3]{ 0.5394782755183589, 0.7809840637246965, 0.3146856883491796 };
	vtkMath::Normalize(adirection);
	sharpplane->SetNormal(adirection);
	auto seamly = Seam(surface, sharpplane);

	std::array<double, 3> tongue_direction {0.0, 1.0, 0.0};
	auto cross_sectly = CrossSection(sharpplane, seamly.get_disks(), tongue_direction);
	auto cross_section_projected = cross_sectly.get_plane_disks();
	write_it(cross_section_projected, "cross_sect_projected.vtp");
	PlanarNormalFilter planer = PlanarNormalFilter(cross_section_projected);
	auto normals = planer.get_planar_normals();
	auto boundwithnormals = planer.get_boundary();
	boundwithnormals->GetPointData()->AddArray(normals);
	boundwithnormals->GetPointData()->AddArray(planer.get_raw_normals());
	write_it(boundwithnormals, "check_dez_normals2.vtp");

	auto offler = planer.offsetter_field(-2.0, 5.0);
	write_it(offler, "hipofffield.vtp");
	auto expandilizer = Expandilizer(offler);
	SeamParameters params;
	auto widget = expandilizer.tongue_and_groove(params);

	vtkSmartPointer<vtkTransform> transform = cross_sectly.get_transform();
	auto inv_transform = vtkSmartPointer<vtkTransform>::New();
	auto inv_transform_matrix = vtkSmartPointer<vtkMatrix4x4>::New();
	transform->GetInverse(inv_transform_matrix);
	inv_transform->SetMatrix(inv_transform_matrix);
	auto inver = vtkSmartPointer<vtkTransformPolyDataFilter>::New();
	inver->SetInputData(widget);
	inver->SetTransform(inv_transform);
	inver->Update();
	auto theseam = inver->GetOutput();
	write_it(theseam, "seam_C.vtp");
}

struct PipelineParameters {
	double groove_outer = 1.0;
	double groove_inner = -1.0;
	double gap_radial = 0.5;
	double gap_depth = 0.5;
	double tongue_depth = 5.0;
	double trim_depth = 1.0;
	std::string input_filepath = "C:\\Users\\sscott\\Pictures\\lattice1.stl";
	std::string name1 = "adam";
};

void pipeline_explorer(PipelineParameters params) {
	/*boost::filesystem::path infilepath{ "C:\\Users\\sscott\\Pictures\\lattice1.stl" };
	auto reader = vtkSmartPointer<vtkSTLReader>::New();
	reader->SetFileName(infilepath.string().c_str());
	reader->Update();
	auto surface = reader->GetOutput();

	auto sharpplane = vtkSmartPointer<vtkPlane>::New();
	double anorigin[3]{ 684.4189787288672, 0.0, 0.0 };
	sharpplane->SetOrigin(anorigin);
	double adirection[3]{ 0.5394782755183589, 0.7809840637246965, 0.3146856883491796 };
	vtkMath::Normalize(adirection);
	sharpplane->SetNormal(adirection);
	auto seamly = Seam(surface, sharpplane);

	std::array<double, 3> tongue_direction{ 0.0, 1.0, 0.0 };
	auto cross_sectly = CrossSection(sharpplane, seamly.get_disks(), tongue_direction);
	auto cross_section_projected = cross_sectly.get_plane_disks();
	write_it(cross_section_projected, "cross_sect_projected.vtp");
	PlanarNormalFilter planer = PlanarNormalFilter(cross_section_projected);
	auto normals = planer.get_planar_normals();
	auto boundwithnormals = planer.get_boundary();
	boundwithnormals->GetPointData()->AddArray(normals);
	boundwithnormals->GetPointData()->AddArray(planer.get_raw_normals());
	write_it(boundwithnormals, "check_dez_normals2.vtp");

	auto offler = planer.offsetter_field(-2.0, 5.0);
	write_it(offler, "hipofffield.vtp");
	auto expandilizer = Expandilizer(offler);
	SeamParameters params;
	auto widget = expandilizer.tongue_and_groove(params);

	vtkSmartPointer<vtkTransform> transform = cross_sectly.get_transform();
	auto inv_transform = vtkSmartPointer<vtkTransform>::New();
	auto inv_transform_matrix = vtkSmartPointer<vtkMatrix4x4>::New();
	transform->GetInverse(inv_transform_matrix);
	inv_transform->SetMatrix(inv_transform_matrix);
	auto inver = vtkSmartPointer<vtkTransformPolyDataFilter>::New();
	inver->SetInputData(widget);
	inver->SetTransform(inv_transform);
	inver->Update();
	auto theseam = inver->GetOutput();
	write_it(theseam, "seam_C.vtp");*/
}



int main() {
	example_lattice_join_different_slide_direction();
};