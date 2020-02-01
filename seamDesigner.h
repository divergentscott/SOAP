#ifndef SEAM_DESIGNER_H
#define SEAM_DESIGNER_H

#ifndef _USE_MATH_DEFINES
#define _USE_MATH_DEFINES
#endif // !1

#include "meshData.h"
#include "planarMeshing.h"

#include <vtkPolyData.h>
#include <vtkSmartPointer.h>

namespace d3d {
	namespace soap {
		namespace point {
			// Porting over the points from the planar namespace in planarMeshing
			using r2 = planar::point::r2;
			using r3 = planar::point::r3;
		};

		class Designer {
		public:
			// Publicly available variables
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
			planar::CrossSectioner::ptr cross_sectioner_;
			planar::PlanarMesh::ptr offset_field_;
			// Above the plane lies the top
			bool is_top_computed;
			vtkSmartPointer<vtkPolyData> top_tongueless_;
			// Below the plane lies the bottom
			bool is_bottom_computed;
			vtkSmartPointer<vtkPolyData> bottom_grooveless_;
			//
			bool is_seam_computed;
			vtkSmartPointer<vtkPolyData> tongue_;
			vtkSmartPointer<vtkPolyData> groove_;

		public:
			//Public methods
			Designer();
			Designer(d3d::CommonMeshData &meshin, Designer::Parameters &parameters);
			Designer(vtkSmartPointer<vtkPolyData> meshin, Designer::Parameters &parameters);
			// Return the groove mesh only
			vtkSmartPointer<vtkPolyData> get_groove();
			// Return the tongue mesh only
			vtkSmartPointer<vtkPolyData> get_tongue();
			// Return the tongue and groove mesh.
			vtkSmartPointer<vtkPolyData> get_seam();
			// Return the top with the tongue geometry
			vtkSmartPointer<vtkPolyData> get_top();
			vtkSmartPointer<vtkPolyData> get_top_tongueless();
			// Return the bottom with the groove geometry
			vtkSmartPointer<vtkPolyData> get_bottom();
			vtkSmartPointer<vtkPolyData> get_bottom_grooveless();
			// Return full jointed geometry
			vtkSmartPointer<vtkPolyData> get_mesh();

		private:
			//Private methods
			void compute_top_tongueless();
			void compute_groove();
			void compute_tongue();
			vtkSmartPointer<vtkPolyData> compute_plane_clip();

		};


	} //end namespace seam
}; //end namespace d3d
#endif