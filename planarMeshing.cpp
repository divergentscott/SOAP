#ifndef _USE_MATH_DEFINES
#define _USE_MATH_DEFINES
#endif // !1

#include "meshData.h"

#include "planarMeshing.h"

#include <vtkSmartPointer.h>
#include <vtkPolyData.h>

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

	CrossSectioner::CrossSectioner() {};

	// From CommonMesh
	CrossSectioner::CrossSectioner(CommonMeshData & mesh, const point::r3 plane_origin, const point::r3 plane_normal) {
		mesh_in_ = makePolydata(mesh);
		CrossSectioner(mesh_in_, plane_origin, plane_normal);
	};

	// From CommonMesh with distinct tongue
	CrossSectioner::CrossSectioner(CommonMeshData & mesh, const point::r3 plane_origin, const point::r3 plane_normal, const point::r3 tongue_direction) {
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


	//class CurveCollection {
	//	//Collection of closed curves in the plane
	//private:
	//	vtkSmartPointer<vtkPolyData> polydata_ = vtkSmartPointer<vtkPolyData>::New();
	//	//The full list of point ids ordered so that each component is oriented and next_point_[foo] follows foo in the oriented traversal of the component. 
	//	std::vector<int> next_point_ = {}
	//	//One point id from each connected component. The number of connected components is the size of the vector.
	//	std::vector<int> component_basepoints_ = {};
	//	// ~Tangent vectors computed at each point. Maybe consider bezier tangents or something.
	//	std::vector<point::r2> tangents_ = {};
	//	// Outward pointing normal vectors. Just tangent vectors rotated 90 degrees and 
	//	std::vector<point::r3> normals_ = {};
	//};

	//class PlanarMesher {
	//private:
	//	vtkSmartPointer<vtkPolyData> polydata_ = vtkSmartPointer<vtkPolyData>::New();

	//};
}; //end namespace planar
}; //end namespace d3d