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
			bool is_planar_transform_computed_ = false;
			vtkSmartPointer<vtkTransform> planar_transform_; //{ 1.0, 0.0, 0.0, 0.0, /**/ 0.0, 1.0, 0.0, 0.0, /**/ 0.0, 0.0, 1.0, 0.0,  /**/ 0.0, 0.0, 0.0, 1.0  /**/ };
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
			// Get the cross section curves
			vtkSmartPointer<vtkPolyData> get_cross_section(const double margin = 0.0);
			// Return the capped surface above or below the plane
			vtkSmartPointer<vtkPolyData> get_clipping(const double margin = 0.0, const bool is_clipping_above = true);
			// Get the transformed cross section curves
			vtkSmartPointer<vtkPolyData> get_planed_cross_section(const double margin = 0.0);
			void compute_planar_transform();
			// Transformation of cut plane into xy plane with sheer to flatten tongue tilt.
			vtkSmartPointer<vtkTransform> get_planar_transform();
		};

		/*~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~*/
		/* Joint Tongue and Groove Geometry Designer */
		/*~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~*/
		// Iterated application needs work!
		// The offsets will interfere with the previous groove.
		// Remeshing is also necessary.

		class SeamDesigner {
		public:
			using ptr = std::shared_ptr<SeamDesigner>;
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
			SeamDesigner::Parameters parameters_;
			vtkSmartPointer<vtkPolyData> vtk_polydata_ = vtkSmartPointer<vtkPolyData>::New();
			vtkSmartPointer<vtkPolyData> offgrid_ = vtkSmartPointer<vtkPolyData>::New();
			vtkSmartPointer<vtkContourFilter> isoer_ = vtkSmartPointer<vtkContourFilter>::New();
			CrossSectioner::ptr cross_sectioner_;
			vtkSmartPointer<vtkPolyData> offset_field_ = vtkSmartPointer<vtkPolyData>::New();
			// Above the plane lies the top
			bool is_top_computed = false;
			vtkSmartPointer<vtkPolyData> top_tongueless_ = vtkSmartPointer<vtkPolyData>::New();
			// Below the plane lies the bottom
			bool is_bottom_computed = false;
			vtkSmartPointer<vtkPolyData> bottom_grooveless_ = vtkSmartPointer<vtkPolyData>::New();
			bool is_seam_computed = false;
			void compute_tongue();
			vtkSmartPointer<vtkPolyData> tongue_ = vtkSmartPointer<vtkPolyData>::New();
			void compute_groove();
			vtkSmartPointer<vtkPolyData> groove_ = vtkSmartPointer<vtkPolyData>::New();
			// Create flat regions from distance field
			vtkSmartPointer<vtkPolyData> flat(const double off, const double z);
			// Create flat regions from distance field
			vtkSmartPointer<vtkPolyData> flat(const double off1, const double off2, const double z);
			// Extrude band from offset curve
			vtkSmartPointer<vtkPolyData> band(const double off, const double z1, const double z2);

		public:
			SeamDesigner();

			SeamDesigner(d3d::CommonMeshData &meshin, SeamDesigner::Parameters &parameters);
			
			SeamDesigner(vtkSmartPointer<vtkPolyData> meshin, SeamDesigner::Parameters &parameters);
			// Compute the top with tongue
			vtkSmartPointer<vtkPolyData> get_top();
			// Compute the bottom with groove
			vtkSmartPointer<vtkPolyData> get_bottom();
			// Return the tongue and groove mesh.
			vtkSmartPointer<vtkPolyData> get_seam();
		};
	} //end namespace seam
}; //end namespace d3d
#endif