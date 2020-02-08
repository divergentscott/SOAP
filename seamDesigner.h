#ifndef SEAM_DESIGNER_H
#define SEAM_DESIGNER_H

#ifndef _USE_MATH_DEFINES
#define _USE_MATH_DEFINES
#endif // !1

#include "meshData.h"
#include "planarMeshing.h"

#include <vtkContourFilter.h>
#include <vtkPolyData.h>
#include <vtkSmartPointer.h>
#include <vtkTransform.h>

namespace d3d {
	namespace soap {
		namespace point {
			using namespace planar::point;
		}

		/*~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~*/
		/* Mesh/Plane Cross Section and Mesh Sides */
		/*~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~*/
		
		class CrossSectioner {
			// This class is intended to slice a cross section from a 2-manifold 
			// Use to  generate cross section curves and a linear transformation 
		public:
			using ptr = std::shared_ptr<CrossSectioner>;
		private:
			vtkSmartPointer<vtkPolyData> mesh_in_ = vtkSmartPointer<vtkPolyData>::New();
			point::r3 plane_origin_{ 0.0, 0.0, 0.0 };
			point::r3 plane_normal_{ 0.0, 0.0, 1.0 };
			bool is_tongue_ = false;
			point::r3 tongue_direction_{ 0.0, 0.0, 0.0 };
			vtkSmartPointer<vtkTransform> planar_transform_;
			bool is_cross_section_computed_ = false;
			vtkSmartPointer<vtkPolyData> cross_section_ = vtkSmartPointer<vtkPolyData>::New();
			vtkSmartPointer<vtkTransform> get_sheer_plane_transform(const point::r3 );

		public:
			CrossSectioner();
			CrossSectioner(const d3d::CommonMeshData & mesh, const point::r3 plane_origin, const point::r3 plane_normal);
			CrossSectioner(vtkSmartPointer<vtkPolyData> mesh, const point::r3 plane_origin, const point::r3 plane_normal);
			// Introducing a tongue_direction will shear the cross section into a plane with the same origin as the cut plane,
			// but the tongue_direction as the normal.
			CrossSectioner(const d3d::CommonMeshData & mesh, const point::r3 plane_origin, const point::r3 plane_normal, const point::r3 tongue_direction);
			CrossSectioner(vtkSmartPointer<vtkPolyData> mesh, const point::r3 plane_origin, const point::r3 plane_normal, const point::r3 tongue_direction);
			// Return the capped surface above or below the plane
			vtkSmartPointer<vtkPolyData> get_clipping(const double margin = 0.0, const bool is_clipping_above = true);
			// Get the cross section curves
			vtkSmartPointer<vtkPolyData> get_cross_section(const double margin = 0.0);
			// Transformation of cut plane into xy plane with sheer to flatten tongue tilt.
			vtkSmartPointer<vtkTransform> get_planar_transform();
		};

		/*~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~*/
		/* Joint Tongue and Groove Geometry Designer */
		/*~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~*/

		class Designer {
		public:
			using ptr = std::shared_ptr<Designer>;
			struct Parameters {
				double groove_outer = 1.0;
				double groove_inner = -1.0;
				double gap_radial = 0.5;
				double gap_depth = 0.5;
				double tongue_depth = 3.0;
				double trim_depth = 1.0;
				point::r3 plane_origin{ 0.0, 0.0, 0.0 };
				point::r3 plane_normal{ 0.0, 0.0, 1.0 };
				bool use_tongue = false;
				point::r3 tongue_direction{ 0.0, 0.0, 0.0 };
			};

		private:
			// Private variables
			Designer::Parameters parameters_;
			vtkSmartPointer<vtkPolyData> vtk_polydata_ = vtkSmartPointer<vtkPolyData>::New();
			vtkSmartPointer<vtkPolyData> offgrid_ = vtkSmartPointer<vtkPolyData>::New();
			vtkSmartPointer<vtkContourFilter> isoer_ = vtkSmartPointer<vtkContourFilter>::New();
			CrossSectioner::ptr cross_sectioner_;
			vtkSmartPointer<vtkPolyData> offset_field_ = vtkSmartPointer<vtkContourFilter>::New();
			// Above the plane lies the top
			bool is_top_computed = false;
			vtkSmartPointer<vtkPolyData> top_tongueless_ = vtkSmartPointer<vtkContourFilter>::New();
			// Below the plane lies the bottom
			bool is_bottom_computed = false;
			vtkSmartPointer<vtkPolyData> bottom_grooveless_ = vtkSmartPointer<vtkContourFilter>::New();
			//
			bool is_seam_computed = false;
			vtkSmartPointer<vtkPolyData> tongue_ = vtkSmartPointer<vtkContourFilter>::New();
			vtkSmartPointer<vtkPolyData> groove_ = vtkSmartPointer<vtkContourFilter>::New();
			//Private methods
			//void compute_top_tongueless();
			//void compute_groove();
			//void compute_tongue();
			//vtkSmartPointer<vtkPolyData> compute_plane_clip();

		public:
			Designer();

			Designer(d3d::CommonMeshData &meshin, Designer::Parameters &parameters);
			
			Designer(vtkSmartPointer<vtkPolyData> meshin, Designer::Parameters &parameters);
			
			// Return the groove mesh only
			//vtkSmartPointer<vtkPolyData> get_groove();
			
			// Return the tongue mesh only
			//vtkSmartPointer<vtkPolyData> get_tongue();
			
			// Return the tongue and groove mesh.
			//vtkSmartPointer<vtkPolyData> get_seam();
			
			// Return the top with the tongue geometry
			//vtkSmartPointer<vtkPolyData> get_top();
			//vtkSmartPointer<vtkPolyData> get_top_tongueless();
			
			// Return the bottom with the groove geometry
			//vtkSmartPointer<vtkPolyData> get_bottom();
			//vtkSmartPointer<vtkPolyData> get_bottom_grooveless();
			// Return full jointed geometry
			//vtkSmartPointer<vtkPolyData> get_mesh();
		};
	} //end namespace seam
}; //end namespace d3d
#endif