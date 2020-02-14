#ifndef PLANAR_MESHING_H
#define PLANAR_MESHING_H

#include <array>
#include <iostream>

#include <vtkSmartPointer.h>
#include <vtkPolyData.h>

#include "meshData.h"

namespace planar {
	namespace point {
		const double epsilon = std::sqrt(std::numeric_limits<double>::epsilon());
		// Points in the plane R^2
		using r2 = std::array<double, 2>;
		r2 operator+(const r2 a, const r2 b);
		r2 operator-(const r2 a, const r2 b);
		r2 operator*(const double a, const r2 b);
		std::ostream & operator << (std::ostream &out, const r2 p);
		double norm(const r2 p);
		r2 normalize(const r2 p);
		double dot(const r2 p, const r2 q);


		// Points in 3-space R^3
		using r3 = std::array<double, 3>;
		r3 operator+(const r3 a, const r3 b);
		r3 operator-(const r3 a, const r3 b);
		r3 operator*(const double a, const r3 b);
		std::ostream & operator << (std::ostream &out, const r3 p);
		double norm(const r3 p);
		r3 normalize(const r3 p);
		double dot(const r3 p, const r3 q);
		r3 cross(const r3 p, const r3 q);
	};

	class CurveCollection {
	public:
		using ptr = std::shared_ptr<CurveCollection>;
	private:
		//Input must represent a collection of disjoint simple closed curves in the plane.
		vtkSmartPointer<vtkPolyData> polydata_;
		bool is_valid_curve_ = true;
		std::vector<point::r2> points_{};
		std::vector<std::array<int, 2>> edges_{};
		// Next point in the traversal
		std::vector<int> order_next_{};
		// Previous point in the traversal
		std::vector<int> order_prev_{};
		// One basepoint for each connected component
		std::vector<int> basepoints_{};
		// Tangent vectors by discrete difference
		std::vector<point::r2> tangents_{};
		void compute_tangents();
		// Normal vectors by orthogonality to tangents
		std::vector<point::r2> normals_{};
		bool is_normals_computed_ = false;
		void compute_normals();
		// Compute a traversal
		bool orient_curves();
		std::array<int, 2> neighborhood_orientation(int pointid);


	public:
		CurveCollection();
		// Poly data must have lines forming closed curves.
		CurveCollection(vtkSmartPointer<vtkPolyData> poly_curve);

		// Retrieve a vtk polydata representation of the curve collection
		void update_polydata();
		vtkSmartPointer<vtkPolyData> get_polydata();

		// Ray to segment intersection computation for ray trace
		point::r2 ray_segment_intersect(const point::r2 orig, const point::r2 dir, const point::r2 a, const point::r2 b);
		bool is_ray_segment_intersect(const point::r2 orig, const point::r2 dir, const point::r2 a, const point::r2 b);
		double distance_till_impact(const int pointid, const point::r2 dir);


		// Compute if a point lies in the interior of the curve collection
		bool is_point_in(const point::r2 point);

		// Write to file for debugging. 
		void write_to_vtp(const std::string filename);

		// Point access
		int get_number_of_points();
		point::r2 get_point(int p_id);

		// Normal vector access
		point::r2 get_normal(int p_id);

		// Signed distance field locally around the curve collection
		vtkSmartPointer<vtkPolyData> distance_field(const double min_dist, const double max_dist, const int sampling=13);

		//todo: edge subdivider, close point merger, signed distance offsetter
	};

}; // end namespace planar
#endif // 