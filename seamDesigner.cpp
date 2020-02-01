#ifndef _USE_MATH_DEFINES
#define _USE_MATH_DEFINES
#endif // !1

#include "seamDesigner.h"

namespace d3d {
	namespace soap {
		Designer::Designer() = default;

		Designer::Designer(d3d::CommonMeshData & meshin, Designer::Parameters & parameters)
		{
			parameters_ = parameters;
			vtk_polydata_ = makePolydata(meshin);
			Designer::Designer(meshin, parameters);
		};

		Designer::Designer(vtkSmartPointer<vtkPolyData> meshin, Designer::Parameters & parameters) {
			parameters_ = parameters;
			vtk_polydata_ = meshin;
			if (parameters_.use_tongue) {
				cross_sectioner_ = std::make_shared<planar::CrossSectioner>(vtk_polydata_,
					parameters_.plane_origin,
					parameters_.plane_normal,
					parameters_.tongue_direction);
			}
			else
			{
				cross_sectioner_ = std::make_shared<planar::CrossSectioner>(vtk_polydata_,
					parameters_.plane_origin,
					parameters_.plane_normal);
			}
			// Compute the cross section.
			// Clean up the curves?
			// Make the seam geometry
			// 
				// Get the offset field
				// Construct the bands and flats
				// Merge all into one mesh
				// Transform seam into place

			// Cut out and cap top
			// Cut out and cap bottom
			// Merge top and tongue.
			// Merge bottom and groove.
			
			// !!!! How do intersection curves fit together?? !!!!
			// Multiple planes needs a method to pass the split sides correctly capped.
		};

		//	auto sharpplane = vtkSmartPointer<vtkPlane>::New();
		//	sharpplane->SetOrigin(params.cut_plane_origin.data());
		//	double adirection[3]{ params.cut_plane_dir[0], params.cut_plane_dir[1], params.cut_plane_dir[2] };
		//	vtkMath::Normalize(adirection);
		//	sharpplane->SetNormal(adirection);
		//	auto seamly = Seam(start_surface, sharpplane);
		//	CrossSection cross_sectly;
		//	if (params.use_tounge) cross_sectly = CrossSection(sharpplane, seamly.get_disks(), params.tongue_dir);
		//	else cross_sectly = CrossSection(sharpplane, seamly.get_disks());
		//	//
		//	auto cross_section_projected = cross_sectly.get_plane_disks();
		//
		//	auto cleaner = vtkSmartPointer<vtkCleanPolyData>::New();
		//	cleaner->PointMergingOn();
		//	cleaner->ConvertLinesToPointsOn();
		//	cleaner->ConvertPolysToLinesOff();
		//	cleaner->SetInputData(cross_section_projected);
		//	//cleaner->SetTolerance(0.0005);
		//	cleaner->Update();
		//	auto clean_disks = cleaner->GetOutput();
		//	write_it(clean_disks, "vtkcleaned.vtp");
		//
		//	PlanarNormalFilter planer = PlanarNormalFilter(clean_disks);
		//	auto normals = planer.get_planar_normals();
		//	auto boundwithnormals = planer.get_boundary();
		//	boundwithnormals->GetPointData()->AddArray(normals);
		//	boundwithnormals->GetPointData()->AddArray(planer.get_raw_normals());
		//	write_it(boundwithnormals, "boundary_normals.vtp");
		//
		//
		//	double z_tongue_tip = (params.gap_depth - params.tongue_depth) / 2.0;
		//	double z_teeth_back = z_tongue_tip - params.gap_depth;
		//
		//	auto offler = planer.offsetter_field((params.groove_inner - params.gap_radial) * 3, params.groove_outer * 3);
		//	write_it(offler, "offler.vtp");
		//	auto expandilizer = Expandilizer(offler);
		//	SeamParameters expan_param;
		//	expan_param.groove_outer = params.groove_outer;
		//	expan_param.groove_inner = params.groove_inner;
		//	expan_param.gap_radial = params.gap_radial;
		//	expan_param.gap_depth = params.gap_depth;
		//	expan_param.tongue_depth = params.tongue_depth;
		//	expan_param.trim_depth = params.trim_depth;
		//
		//	auto widget = expandilizer.tongue_and_groove(expan_param);
		//
		//	vtkSmartPointer<vtkTransform> transform = cross_sectly.get_transform();
		//	auto inv_transform = vtkSmartPointer<vtkTransform>::New();
		//	auto inv_transform_matrix = vtkSmartPointer<vtkMatrix4x4>::New();
		//	transform->GetInverse(inv_transform_matrix);
		//	inv_transform->SetMatrix(inv_transform_matrix);
		//	auto inver = vtkSmartPointer<vtkTransformPolyDataFilter>::New();
		//	inver->SetInputData(widget);
		//	inver->SetTransform(inv_transform);
		//	inver->Update();
		//	auto theseam = inver->GetOutput();
		//	auto seam_polydatas = get_components(theseam);
		//	int cnt = 2;
		//	for (auto pol : seam_polydatas) {
		//		write_it(pol, params.name1 + std::to_string(cnt) + ".vtp");
		//		cnt++;
		//	}
		//	auto clipper = vtkSmartPointer<vtkClipPolyData>::New();
		//	clipper->SetClipFunction(sharpplane);
		//	clipper->SetInputData(start_surface);
		//
		//	clipper->SetValue(-z_teeth_back);
		//	clipper->Update();
		//	write_it(clipper->GetOutput(), params.name1 + "0.vtp");
		//	clipper->SetValue(z_teeth_back);
		//	clipper->GenerateClippedOutputOn();
		//	clipper->Update();
		//	write_it(clipper->GetClippedOutput(), params.name1 + "1.vtp");
		//}
	}; // end namespace seam
}; //end namespace d3d