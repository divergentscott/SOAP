#ifndef PLANAR_MESHING_H
#define PLANAR_MESHING_H

#include <array>
#include <iostream>

#include <vtkSmartPointer.h>
#include <vtkPolyData.h>
#include "meshData.h"

namespace d3d {
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
			std::vector<point::r2> points_ {};
			// Next point in the traversal
			std::vector<int> point_next_ {};
			// Previous point in the traversal
			std::vector<int> point_prev_ {};
			// One basepoint for each connected component
			std::vector<int> basepoints_{};
			std::vector<point::r2> tangents_ {};
			std::vector<point::r2> normals_ {};
			vtkSmartPointer<vtkPolyData> poly_curve_;
			bool orient_curves();
			void compute_tangents();
			void compute_normals();

		public:
			CurveCollection();
			//Poly data must have lines forming closed curves.
			CurveCollection(vtkSmartPointer<vtkPolyData> poly_curve);
			point::r2 ray_segment_intersect(point::r2 orig, point::r2 r, point::r2 a, point::r2 b);
			bool is_ray_segment_intersect(point::r2 orig, point::r2 r, point::r2 a, point::r2 b);
			bool is_point_in(point::r2 point);
		};

		class PlanarMesh {
		public:
			using ptr = std::shared_ptr<PlanarMesh>;
		private:
			double dip;
		public:
			PlanarMesh() {};
			//PlanarMesh(vtkSmartPointer<vtkPolyData> mesh);
			//remesh();
			//CurveCollection::ptr get_boundary();
			//set_scalar_field();?
			//get_scalar_field();?
			//CurveCollection::ptr get_isocurves();
		};

		class CrossSectioner {
			// This class is intended to slice a cross section from a 2-manifold 
			// Use to  generate cross section curves and a linear transformation 
		public:
			using ptr = std::shared_ptr<CrossSectioner>;
		private:
			vtkSmartPointer<vtkPolyData> mesh_in_;
			point::r3 plane_origin_{ 0.0, 0.0, 0.0 };
			point::r3 plane_normal_ { 0.0, 0.0, 1.0 };
			point::r3 tongue_direction_{ 0.0, 0.0, 0.0 };

		public:
			CrossSectioner();
			CrossSectioner(d3d::CommonMeshData & mesh, const point::r3 plane_origin, const point::r3 plane_normal);
			CrossSectioner(vtkSmartPointer<vtkPolyData> mesh, const point::r3 plane_origin, const point::r3 plane_normal);
			// Introducing a tongue_direction will shear the cross section into a plane with the same origin as the cut plane,
			// but the tongue_direction as the normal.
			CrossSectioner(d3d::CommonMeshData & mesh, const point::r3 plane_origin, const point::r3 plane_normal, const point::r3 tongue_direction);
			CrossSectioner(vtkSmartPointer<vtkPolyData> mesh, const point::r3 plane_origin, const point::r3 plane_normal, const point::r3 tongue_direction);
			// Return the capped surface above or below the plane
			vtkSmartPointer<vtkPolyData> get_clipping(double margin=0.0, bool is_clipping_above = true);
			vtkSmartPointer<vtkPolyData> get_cross_section(double margin = 0.0);
		};

}; // end namespace planar
}; //end namespace d3d
#endif // !1