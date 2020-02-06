#include "meshData.h"

#include "planarMeshing.h"

#include <vtkAppendPolyData.h>
#include <vtkCleanPolyData.h>
#include <vtkContourTriangulator.h>
#include <vtkClipPolyData.h>
#include <vtkCutter.h>
#include <vtkDoubleArray.h>
#include <vtkIdList.h>
#include <vtkLine.h>
#include <vtkPlane.h>
#include <vtkPointData.h>
#include <vtkPolyData.h>
#include <vtkSmartPointer.h>
#include <vtkXMLPolyDataWriter.h>

namespace d3d{
namespace planar {

	namespace point {
		// Points in the plane R^2
		r2 operator+(const r2 a, const r2 b) {
			r2 c{ a[0] + b[0], a[1] + b[1] };
			return c;
		}

		r2 operator-(const r2 a, const r2 b) {
			r2 c{ a[0] - b[0], a[1] - b[1] };
			return c;
		}

		r2 operator*(const double a, const r2 b) {
			r2 c{ a*b[0], a*b[1] };
			return c;
		}

		ostream & operator << (ostream &out, const r2 p) {
			out << p[0] << " " << p[1];
			return out;
		}

		double norm(const r2 p) {
			return std::sqrt(p[0] * p[0] + p[1] * p[1]);
		}

		r2 normalize(const r2 p) {
			return (1 / norm(p)) * p;
		}

		double dot(const r2 p, const r2 q) {
			return p[0] * q[0] + p[1] * q[1];
		};
		//end R^2 point operations 

		// Points in 3-space R^3
		r3 operator+(const r3 a, const r3 b) {
			r3 c{ a[0] + b[0], a[1] + b[1], a[2] + b[2] };
			return c;
		}

		r3 operator-(const r3 a, const r3 b) {
			r3 c{ a[0] - b[0], a[1] - b[1], a[2] - b[2] };
			return c;
		}

		r3 operator*(const double a, const r3 b) {
			r3 c{ a*b[0], a*b[1] , a*b[2] };
			return c;
		}

		ostream & operator << (ostream &out, const r3 p) {
			out << p[0] << " " << p[1] << " " << p[2];
			return out;
		}

		double norm(const r3 p) {
			return std::sqrt(p[0] * p[0] + p[1] * p[1] + p[2] * p[2]);
		}

		r3 normalize(const r3 p) {
			return (1 / norm(p)) * p;
		}

		double dot(const r3 p, const r3 q) {
			return p[0] * q[0] + p[1] * q[1] + p[2] * q[2];
		};

		r3 cross(const r3 p, const r3 q) {
			return { p[1] * q[2] - q[2] * p[1], p[2] * q[0] - q[0] * p[2] , p[0] * q[1] - q[1] * p[0] };
		};

		//end R^3 point operations
		template<class T>
		bool within_epsilon(T p, T q) {
			return (norm(p - q) < epsilon);
		}

	}; //end namespace point
	using namespace point;

	// Tao is 2pi
	const double CONSTANT_TAO = 6.28318530717958647692528676655900576839433879875021164194;

	// Random number generator
	double random_angle() {
		/*
		Random double in the range (-0.5,0.5)
		*/
		return ((double)std::rand()) / RAND_MAX * CONSTANT_TAO;
	}

	// Cross Sectioner
	CrossSectioner::CrossSectioner() {};

	// From CommonMesh
	CrossSectioner::CrossSectioner(const CommonMeshData & mesh, const point::r3 plane_origin, const point::r3 plane_normal) {
		mesh_in_ = makePolydata(mesh);
		CrossSectioner(mesh_in_, plane_origin, plane_normal);
	};

	// From CommonMesh with distinct tongue
	CrossSectioner::CrossSectioner(const CommonMeshData & mesh, const point::r3 plane_origin, const point::r3 plane_normal, const point::r3 tongue_direction) {
		mesh_in_ = makePolydata(mesh);
		CrossSectioner(mesh_in_, plane_origin, plane_normal, tongue_direction);
	};

	// From Polydata
	CrossSectioner::CrossSectioner(vtkSmartPointer<vtkPolyData> mesh, const point::r3 plane_origin, const point::r3 plane_normal) {
		mesh_in_ = mesh;
	};

	// From Polydata with distinct tongue
	CrossSectioner::CrossSectioner(vtkSmartPointer<vtkPolyData> mesh, const point::r3 plane_origin, const point::r3 plane_normal, const point::r3 tongue_direction) {
		//plane_normal_ = point::normalize(plane_normal);
		double tongue_norm = point::norm(tongue_direction);
		// If the tongue_direction is too small or too close to the plane_normal, ignore it. Otherwise the sheer transformation is poorly conditioned.
		if (tongue_norm < point::epsilon) {
			CrossSectioner::CrossSectioner(mesh, plane_origin, plane_normal);
		}
		else {
			tongue_direction_ = (1.0/tongue_norm) * tongue_direction;
			mesh_in_ = mesh;
		}
	};

	// Return the capped surface above the plane
	vtkSmartPointer<vtkPolyData> CrossSectioner::get_clipping(const double margin, const bool is_clipping_above) {
		//First get the surface above/below the margin
		auto clipper = vtkSmartPointer<vtkClipPolyData>::New();
		clipper->SetInputData(mesh_in_);
		auto cutting_plane = vtkSmartPointer<vtkPlane>::New();
		cutting_plane->SetOrigin(plane_origin_.data());
		cutting_plane->SetNormal(plane_normal_.data());
		clipper->SetClipFunction(cutting_plane);
		clipper->SetInsideOut(!is_clipping_above);
		clipper->SetValue(margin);
		clipper->Update();

		// Get the cross section at the margin and fill it to cap the surface
		// Cut cross-section curves
		auto cutter = vtkSmartPointer<vtkCutter>::New();
		cutter->SetInputData(mesh_in_);
		cutter->SetCutFunction(cutting_plane);
		cutter->SetValue(0, margin);
		cutter->Update();

		// triangulate cross-section curves.
		auto filler = vtkSmartPointer<vtkContourTriangulator>::New();
		filler->SetInputConnection(cutter->GetOutputPort());
		filler->Update();

		//Combine
		auto appender = vtkSmartPointer<vtkAppendPolyData>::New();
		appender->AddInputConnection(clipper->GetOutputPort());
		appender->AddInputConnection(cutter->GetOutputPort());
		appender->Update();

		//Clean
		auto cleaner = vtkSmartPointer<vtkCleanPolyData>::New();
		cleaner->AddInputConnection(appender->GetOutputPort());
		cleaner->PointMergingOn();
		cleaner->Update();

		return cleaner->GetOutput();
	};
	
	vtkSmartPointer<vtkPolyData> CrossSectioner::get_cross_section(const double margin) {
		//The cutting plane
		auto cutting_plane = vtkSmartPointer<vtkPlane>::New();
		cutting_plane->SetOrigin(plane_origin_.data());
		cutting_plane->SetNormal(plane_normal_.data());

		// Cut cross-section curves
		auto cutter = vtkSmartPointer<vtkCutter>::New();
		cutter->SetInputData(mesh_in_);
		cutter->SetCutFunction(cutting_plane);
		cutter->SetValue(0, margin);
		cutter->Update();
		return cutter->GetOutput();
	};

	/* COLLECTION OF SIMPLE CLOSED CURVES IN THE PLANE */
	point::r2 CurveCollection::ray_segment_intersect(const point::r2 orig, const point::r2 dir, const point::r2 a, const point::r2 b) {
		// Compute ray, line segment intersection.
		// Ray has origin orig and direction r. That is: { orig + t * r | t>0 } 
		// Line segment endpoints a b. That is: { a * s + b * (1-s) | s in [0,1] }
		// Computes the parameter t and s at intersection.
		// If the ray and segment are colinear, there may be infinitely many solutions.
		// In that case return the s = 0 solution, so t is such that orig + t * r = b.
		//  orign + t * r = 
		double a01 = b[0] - a[0];
		double a11 = b[1] - a[1];
		double det = dir[0] * a11 - a01 * dir[1];
		if (std::abs(det) > point::epsilon) {
			//Noncolinear
			double s0 = b[0] - orig[0];
			double s1 = b[1] - orig[1];
			double t0 = (a11*s0 - a01 * s1) / det;
			double t1 = (dir[0] * s1 - dir[1] * s0) / det;
			return { t0,t1 };
		}
		else {
			//Colinear
			double t1;
			if (a01 > point::epsilon) {
				t1 = (b[0] - orig[0]) / a01;
			}
			else {
				t1 = (b[1] - orig[1]) / a11;
			}
			if ((-point::epsilon < t1) & (t1 < 1 + point::epsilon)) {
				return { 0.0, t1 };
			}
			else {
				if (dir[0] > point::epsilon) {
					return { (b[0] - orig[0]) / dir[0], 0.0 };
				}
				else {

					return { (b[1] - orig[1]) / dir[1], 0.0 };
				}
			}
		}
	}

	bool CurveCollection::is_ray_segment_intersect(const point::r2 orig, const point::r2 dir, const point::r2 a, const point::r2 b) {
		point::r2 t = ray_segment_intersect(orig, dir, a, b);
		return (t[0] > -point::epsilon) & (-point::epsilon < t[1]) & (t[1] < 1 + point::epsilon);
	}

	bool CurveCollection::is_point_in(point::r2 point) {
		/* Determine if a point lies in a curve collecion via parity of a random direction ray trace. */
		double theta = random_angle();
		point::r2 direct{ cos(theta), sin(theta) };
		int isect_count = 0;
		for (int foo = 0; foo < edges_.size(); foo++) {
			point::r2 point_foo = points_[edges_[foo][0]];
			point::r2 point_next = points_[edges_[foo][1]];
			std::array<double, 2> isect_params = ray_segment_intersect(point, direct, point_foo, point_next);
			// Failure!
			if ((isect_params[1] == 0.0)&(isect_params[0] > 0.0)) {
				// If an edge is a subset of a ray, the pairity test will fail.
				// The only hope to resolve is to check the direction of  nonparallel edges o 
				// This is slightly dangerous! Technically could be an infinite loop...
				// If the segment fully lies inside the ray, we reroll the ray direction and restart.
				// std::cout << "Segment is subset of ray! Reroll";
				theta = random_angle();
				direct = { cos(theta), sin(theta) };
				isect_count = 0;
				foo = -1;
				continue;
			}
			// But normally the ray hits the segment interior.
			if ((-point::epsilon < isect_params[1]) & (isect_params[1] < 1 + point::epsilon)) {
				if (std::abs(isect_params[0]) <= point::epsilon) {
					return true;
				}
				else if (isect_params[0] > 0) {
					isect_count++;
				}
			}
		}
		return (isect_count % 2) == 1;
	}

	double CurveCollection::distance_till_impact(const int pointid, const point::r2 dir) {
		/* Compute the distance on a line at orig with direction dir travels in either direction before hitting an edge not containing pointid.*/
		double dist = std::numeric_limits<double>::infinity();
		point::r2 at = points_[pointid];
		for (auto edge : edges_) {
			if ((pointid != edge[0]) & (pointid != edge[1])) {
				point::r2 impact = ray_segment_intersect(at, dir, points_[edge[0]], points_[edge[1]]);
				//impact[0] is the length along the ray, impact[1] is the length along the line containing edge
				if ((impact[1] < 1 + point::epsilon) & (impact[1] > -point::epsilon)) {
					dist = std::min(std::abs(impact[0]), dist);
				}
			}
		}
		return dist;
	}

	std::array<int, 2> CurveCollection::neighborhood_orientation(int pointid) {
		/*
		 Compute curve orientation at a point by determining if the bordism lies to the left or right.
		 Return <previous point, next point>
		 Ray traces against all edges from at least 3 points.
		*/
		// Construct a point to the left to test if that's inside the curve
		point::r2 at = points_[pointid];
		int p0 = order_next_[pointid];
		int p1 = order_prev_[pointid];
		point::r2 tan = point::normalize(points_[p0] - points_[p1]);
		point::r2 check_dir{ -tan[1], tan[0] };
		// If nearby edges are closer then then check_len, there's no guarantee of correctness!
		double impact_dist = distance_till_impact(pointid, check_dir);
		std::cout << "Impact: " << impact_dist << "\n";
		double check_len = 0.5 * std::min(1.0, impact_dist);
		bool is_left_in, is_right_in;
		do {
			// Ensure that the point is one-sided. 
			point::r2 check_left = at + check_len * check_dir;
			point::r2 check_right = at - check_len * check_dir;
			is_left_in = is_point_in(check_left);
			is_right_in = is_point_in(check_right);
			check_len *= 0.5;
		} while ( (is_left_in==is_right_in) & (check_len > point::epsilon) );	
		if (is_left_in) {
			// If the left side is inward, keep the orientation
			return { p1, p0 };
		}
		else {
			// If the right side is inward, swap the orientation
			return { p0, p1 };
		}
	}

	bool CurveCollection::orient_curves() {
		/*
		Compute the next and previous point in a traversal of all curves.
		*/
		poly_curve_->BuildLinks();
		bool is_valid_curve = true;
		// First initialize the orientation by the edge list.
		order_next_.resize(edges_.size());
		order_prev_.resize(edges_.size());
		for (int foo = 0; foo < order_next_.size(); foo++) {
			order_next_[edges_[foo][0]] = edges_[foo][1];
			order_prev_[edges_[foo][1]] = edges_[foo][0];
		}
		// Traverse the points and swap the incorrectly oriented edges.
		std::vector<bool> point_done(points_.size(), false);
		for (int foo = 0; foo < points_.size(); foo++) {
			if (point_done[foo]) continue;
			else {
				int basepoint = foo;
				basepoints_.push_back(basepoint);
				int p_was = basepoint;
				// Orient so that inside lies on the left of the curve.
				std::array<int, 2> orient = neighborhood_orientation(basepoint);
				order_prev_[basepoint] = orient[0];
				order_next_[basepoint] = orient[1];
				int p_at = order_next_[basepoint];
				point_done[basepoint] = true;
				// Flow around the connected component, correcting the direction.
				while (p_at != basepoint) {
					if (order_prev_[p_at] == p_was) {
						point_done[p_at] = true;
						p_was = p_at;
						p_at = order_next_[p_at];
					}
					else {
						//correct the next point.
						int true_next = order_prev_[p_at];
						order_next_[p_at] = true_next;
						order_prev_[p_at] = p_was;
						point_done[p_at] = true;
						//increment
						p_was = p_at;
						p_at = true_next;
					}
				}
			}
		}
		return is_valid_curve;
	}

	void CurveCollection::compute_tangents() {
		tangents_.resize(points_.size());
		for (int foo = 0; foo < points_.size(); foo++) {
			point::r2 pnext = points_[order_next_[foo]];
			point::r2 pprev = points_[order_prev_[foo]];
			point::r2 tangent_unnormalized = pnext-pprev;
			tangents_[foo] = normalize(tangent_unnormalized);
		}
	}

	void CurveCollection::compute_normals(){
		compute_tangents();
		normals_.resize(points_.size());
		for (int foo = 0; foo < points_.size(); foo++) {
			point::r2 normal{tangents_[foo][1], -tangents_[foo][0]};
			normals_[foo] = normal;
		}
	}

	CurveCollection::CurveCollection(vtkSmartPointer<vtkPolyData> poly_curve)
	{
		poly_curve_ = poly_curve;
		// Size allocation
		points_.resize(poly_curve_->GetNumberOfPoints());
		edges_.resize(poly_curve_->GetNumberOfPoints());
		// Copy the points
		for (int foo = 0; foo < poly_curve_->GetNumberOfPoints(); foo++) {
			double tmp_point_coords[3];
			poly_curve_->GetPoint(foo, tmp_point_coords);
			points_[foo] = { tmp_point_coords[0], tmp_point_coords[1] };
		}
		bool is_valid_curve = true;
		//Copy the edges
		for (int foo = 0; foo < poly_curve_->GetNumberOfLines(); foo++) {
			auto cell_pts = vtkSmartPointer<vtkIdList>::New();
			poly_curve_->GetCellPoints(foo, cell_pts);
			if (cell_pts->GetNumberOfIds() != 2) is_valid_curve = false;
			std::array<int, 2> edge_foo;
			edge_foo[0] = cell_pts->GetId(0);
			edge_foo[1] = cell_pts->GetId(1);
			edges_[foo] = edge_foo;
		}
		bool is_orientable = orient_curves();
		is_valid_curve = is_valid_curve & is_orientable;
	}

	void CurveCollection::write_to_vtp(const std::string outputfilename) {
		/* Debugging output writer for evaluating the point normals. */
		compute_normals();
		auto out_data = vtkSmartPointer<vtkPolyData>::New();
		auto out_points = vtkSmartPointer<vtkPoints>::New();
		auto out_lines = vtkSmartPointer<vtkCellArray>::New();
		auto out_normals = vtkSmartPointer<vtkDoubleArray>::New();
		out_normals->SetName("PlanarNormals");
		out_normals->SetNumberOfComponents(2);
		for (int foo = 0; foo < points_.size(); foo++) {
			double pt_foo[3]{points_[foo][0], points_[foo][1], 0.0};
			out_points->InsertNextPoint(pt_foo);
			auto lien = vtkSmartPointer<vtkLine>::New();
			lien->GetPointIds()->SetId(0, foo);
			lien->GetPointIds()->SetId(1, order_next_[foo]);
			out_lines->InsertNextCell(lien);
			out_normals->InsertNextTuple2(normals_[foo][0], normals_[foo][1]);
		}
		out_data->SetPoints(out_points);
		out_data->SetLines(out_lines);
		out_data->GetPointData()->AddArray(out_normals);
		auto writer = vtkSmartPointer<vtkXMLPolyDataWriter>::New();
		writer->SetFileName(outputfilename.c_str());
		writer->SetInputData(out_data);
		writer->Update();
	}
	//};
}; //end namespace planar
}; //end namespace d3d