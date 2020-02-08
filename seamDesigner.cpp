#ifndef _USE_MATH_DEFINES
#define _USE_MATH_DEFINES
#endif // !1

#include "seamDesigner.h"
#include <vtkAppendPolyData.h>
#include <vtkCenterOfMass.h>
#include <vtkCleanPolyData.h>
#include <vtkClipPolyData.h>
#include <vtkContourTriangulator.h>
#include <vtkCutter.h>
#include <vtkPlane.h>

namespace d3d {
	namespace soap {
		using namespace point;
		/*~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~*/
		/* MESH/PLANE CROSS SECTION AND MESH SIDES */
		/*~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~*/

		// Cross Sectioner
		CrossSectioner::CrossSectioner() = default;

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
			plane_origin_ = plane_origin;
			plane_normal_ = point::normalize(plane_normal);
		};

		// From Polydata with distinct tongue
		CrossSectioner::CrossSectioner(vtkSmartPointer<vtkPolyData> mesh, const point::r3 plane_origin, const point::r3 plane_normal, const point::r3 tongue_direction) {
			//plane_normal_ = point::normalize(plane_normal);
			double tongue_norm = point::norm(tongue_direction);
			// If the tongue_direction is too small or too close to the plane_normal, ignore it. Otherwise the sheer transformation is poorly conditioned.
			if (tongue_norm < point::epsilon) {
				is_tongue_ = false;
				CrossSectioner::CrossSectioner(mesh, plane_origin, plane_normal);
			}
			else {
				is_tongue_ = true;
				tongue_direction_ = (1.0 / tongue_norm) * tongue_direction;
				mesh_in_ = mesh;
				plane_origin_ = plane_origin;
				plane_normal_ = point::normalize(plane_normal);
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

		//
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
		}

		vtkSmartPointer<vtkTransform> CrossSectioner::get_sheer_plane_transform(point::r3 origin) {
			// If the slide direction is not the cut plane direction, then you want to sheer the cross section into
			// a plane whose normal is slide direction, design the joint, the resheer it back in the cut-plane.
			auto smosher = vtkSmartPointer<vtkTransform>::New();
			smosher->PostMultiply();
			if (!is_tongue_) {
				return smosher;
			}
			else {
				point::r3 dircross = point::cross(plane_normal_, tongue_direction_);
				double sin_dir_angle = point::norm(dircross);
				double cos_dir_angle = point::dot(plane_normal_, tongue_direction_);
				dircross = (1.0 / sin_dir_angle) * dircross;
				double dir_basis[3][3]{ {tongue_direction_[0], plane_normal_[0], dircross[0]}, {tongue_direction_[1], plane_normal_[1], dircross[1]}, {tongue_direction_[2], plane_normal_[2], dircross[2]} };
				double dir_basis_inv[3][3];
				vtkMath::Invert3x3(dir_basis, dir_basis_inv);
				double tar_basis[3][3]{ {0.0, sin_dir_angle, 0.0},{0.0,0.0,1.0},{1.0,cos_dir_angle,0.0} };
				double qmat[3][3], qmat_inv[3][3], sheer_mat[3][3], qmat_tmp[3][3];
				double sheer_mat_lower[3][3]{ {1.0, 0.0, 0.0},{0.0, 1.0, 0.0}, {sin_dir_angle / cos_dir_angle, 0.0, 1.0} };
				vtkMath::Multiply3x3(tar_basis, dir_basis_inv, qmat_tmp);
				vtkMath::Orthogonalize3x3(qmat_tmp, qmat);
				vtkMath::Transpose3x3(qmat, qmat_inv);
				double tempmat[3][3];
				vtkMath::Multiply3x3(qmat_inv, sheer_mat_lower, tempmat);
				vtkMath::Multiply3x3(tempmat, qmat, sheer_mat);
				smosher->Translate(-origin[0], -origin[1], -origin[2]);
				double sheer_mat_elements[16]{ sheer_mat[0][0], sheer_mat[0][1], sheer_mat[0][2], 0.0, sheer_mat[1][0], sheer_mat[1][1], sheer_mat[1][2], 0.0, sheer_mat[2][0], sheer_mat[2][1], sheer_mat[2][2], 0.0, 0.0, 0.0, 0.0,1.0 };
				smosher->Concatenate(sheer_mat_elements);
				smosher->Translate(origin[0], origin[1], origin[2]);
				return smosher;
			}
		}

		//
		vtkSmartPointer<vtkTransform> CrossSectioner::get_planar_transform() {
			auto cs = get_cross_section();
			auto transform = vtkSmartPointer<vtkTransform>::New();
			if (cs->GetNumberOfPoints() == 0) {
				return transform;
			}
			else {
				//Compute cut center of mass as new origin
				double origin_tmp0[3];
				vtkCenterOfMass::ComputeCenterOfMass(cs->GetPoints(), nullptr, origin_tmp0);
				point::r3 origin{ origin_tmp0[0], origin_tmp0[1], origin_tmp0[2] };
				origin = origin - point::dot(origin, plane_normal_) * plane_normal_;
				auto transform = get_sheer_plane_transform(origin);
				//auto ppoints = vtkSmartPointer<vtkPoints2D>::New();
				//Project all points of the cross section polys into the plane with the {planarx_, planary_} orthonormal basis
				double tmp_Tpoint[3];
				double tmp_point[3];
				cs->GetPoint(0, tmp_point);
				transform->TransformPoint(tmp_point, tmp_Tpoint);
				double planarx[3], planary[3];
				vtkMath::Subtract(tmp_Tpoint, origin.data(), planarx);
				vtkMath::Normalize(planarx);
				vtkMath::Cross(tongue_direction_.data(), planarx, planary);
				vtkMath::Normalize(planary);
				//The ppoints are the in-plane points: planar points.
				transform->PostMultiply();
				transform->Translate(-origin[0], -origin[1], -origin[2]);
				const double ortho[16]{ planarx[0],planarx[1],planarx[2],0.0, planary[0],planary[1],planary[2],0.0, tongue_direction_[0], tongue_direction_[1], tongue_direction_[2],0.0,  0.0,0.0,0.0,1.0 };
				transform->Concatenate(ortho);
				return transform;
			}
		};
//
//		/*~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~*/
//		/* Joint Tongue and Groove Geometry Designer */
//		/*~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~*/
//
		//Designer::Designer() = default;
//
//		Designer::Designer(d3d::CommonMeshData & meshin, Designer::Parameters & parameters)
//		{
//			parameters_ = parameters;
//			vtk_polydata_ = makePolydata(meshin);
//			Designer(meshin, parameters);
//		};
//
//		Designer::Designer(vtkSmartPointer<vtkPolyData> meshin, Designer::Parameters & parameters) {
//			parameters_ = parameters;
//			vtk_polydata_ = meshin;
//			if (parameters_.use_tongue) {
//				cross_sectioner_ = std::make_shared<CrossSectioner>(vtk_polydata_,
//					parameters_.plane_origin,
//					parameters_.plane_normal,
//					parameters_.tongue_direction);
//			}
//			else
//			{
//				cross_sectioner_ = std::make_shared<CrossSectioner>(vtk_polydata_,
//					parameters_.plane_origin,
//					parameters_.plane_normal);
//			}
//			// Compute the cross section.
//			auto cs = cross_sectioner_->get_cross_section();
//			//planar::CurveCollection::ptr cross_section 
//			// Clean up the curves?
//			// Make the seam geometry
//
//				// Get the offset field
//				// Construct the bands and flats
//				// Merge all into one mesh
//				// Transform seam into place
//
//			// Cut out and cap top
//			// Cut out and cap bottom
//			// Merge top and tongue.
//			// Merge bottom and groove.
//			
//			// !!!! How do intersection curves fit together?? !!!!
//			// Multiple planes needs a method to pass the split sides correctly capped.
//		};
//
//		//	auto sharpplane = vtkSmartPointer<vtkPlane>::New();
//		//	sharpplane->SetOrigin(params.cut_plane_origin.data());
//		//	double adirection[3]{ params.cut_plane_dir[0], params.cut_plane_dir[1], params.cut_plane_dir[2] };
//		//	vtkMath::Normalize(adirection);
//		//	sharpplane->SetNormal(adirection);
//		//	auto seamly = Seam(start_surface, sharpplane);
//		//	CrossSection cross_sectly;
//		//	if (params.use_tounge) cross_sectly = CrossSection(sharpplane, seamly.get_disks(), params.tongue_dir);
//		//	else cross_sectly = CrossSection(sharpplane, seamly.get_disks());
//		//	//
//		//	auto cross_section_projected = cross_sectly.get_plane_disks();
//		//
//		//	auto cleaner = vtkSmartPointer<vtkCleanPolyData>::New();
//		//	cleaner->PointMergingOn();
//		//	cleaner->ConvertLinesToPointsOn();
//		//	cleaner->ConvertPolysToLinesOff();
//		//	cleaner->SetInputData(cross_section_projected);
//		//	//cleaner->SetTolerance(0.0005);
//		//	cleaner->Update();
//		//	auto clean_disks = cleaner->GetOutput();
//		//	write_it(clean_disks, "vtkcleaned.vtp");
//		//
//		//	PlanarNormalFilter planer = PlanarNormalFilter(clean_disks);
//		//	auto normals = planer.get_planar_normals();
//		//	auto boundwithnormals = planer.get_boundary();
//		//	boundwithnormals->GetPointData()->AddArray(normals);
//		//	boundwithnormals->GetPointData()->AddArray(planer.get_raw_normals());
//		//	write_it(boundwithnormals, "boundary_normals.vtp");
//		//
//		//
//		//	double z_tongue_tip = (params.gap_depth - params.tongue_depth) / 2.0;
//		//	double z_teeth_back = z_tongue_tip - params.gap_depth;
//		//
//		//	auto offler = planer.offsetter_field((params.groove_inner - params.gap_radial) * 3, params.groove_outer * 3);
//		//	write_it(offler, "offler.vtp");
//		//	auto expandilizer = Expandilizer(offler);
//		//	SeamParameters expan_param;
//		//	expan_param.groove_outer = params.groove_outer;
//		//	expan_param.groove_inner = params.groove_inner;
//		//	expan_param.gap_radial = params.gap_radial;
//		//	expan_param.gap_depth = params.gap_depth;
//		//	expan_param.tongue_depth = params.tongue_depth;
//		//	expan_param.trim_depth = params.trim_depth;
//		//
//		//	auto widget = expandilizer.tongue_and_groove(expan_param);
//		//
//		//	vtkSmartPointer<vtkTransform> transform = cross_sectly.get_transform();
//		//	auto inv_transform = vtkSmartPointer<vtkTransform>::New();
//		//	auto inv_transform_matrix = vtkSmartPointer<vtkMatrix4x4>::New();
//		//	transform->GetInverse(inv_transform_matrix);
//		//	inv_transform->SetMatrix(inv_transform_matrix);
//		//	auto inver = vtkSmartPointer<vtkTransformPolyDataFilter>::New();
//		//	inver->SetInputData(widget);
//		//	inver->SetTransform(inv_transform);
//		//	inver->Update();
//		//	auto theseam = inver->GetOutput();
//		//	auto seam_polydatas = get_components(theseam);
//		//	int cnt = 2;
//		//	for (auto pol : seam_polydatas) {
//		//		write_it(pol, params.name1 + std::to_string(cnt) + ".vtp");
//		//		cnt++;
//		//	}
//		//	auto clipper = vtkSmartPointer<vtkClipPolyData>::New();
//		//	clipper->SetClipFunction(sharpplane);
//		//	clipper->SetInputData(start_surface);
//		//
//		//	clipper->SetValue(-z_teeth_back);
//		//	clipper->Update();
//		//	write_it(clipper->GetOutput(), params.name1 + "0.vtp");
//		//	clipper->SetValue(z_teeth_back);
//		//	clipper->GenerateClippedOutputOn();
//		//	clipper->Update();
//		//	write_it(clipper->GetClippedOutput(), params.name1 + "1.vtp");
//		//}
//
//		/*
//		Expandilizer(vtkSmartPointer<vtkPolyData> offgrid) {
//			offgrid_ = offgrid;
//			write_it(offgrid_, "offgrid_was.vtp");
//			offgrid_->GetPointData()->SetActiveScalars("offset");
//			isoer_->SetInputData(offgrid_);
//			isoer_->ComputeScalarsOff();
//		}
//
//		vtkSmartPointer<vtkPolyData> flat(double off, double z) {
//			isoer_->SetNumberOfContours(1);
//			isoer_->SetValue(0, off);
//			isoer_->Update();
//			auto filler = vtkSmartPointer<vtkContourTriangulator>::New();
//			filler->SetInputConnection(isoer_->GetOutputPort());
//			filler->Update();
//			auto pprojection = vtkSmartPointer<vtkTransform>::New();
//			pprojection->Translate(0.0, 0.0, z);
//			auto transform_filter = vtkSmartPointer<vtkTransformPolyDataFilter>::New();
//			transform_filter->SetTransform(pprojection);
//			transform_filter->SetInputData(filler->GetOutput());
//			transform_filter->Update();
//			return transform_filter->GetOutput();
//		}
//
//		vtkSmartPointer<vtkPolyData> flat(double off1, double off2, double z) {
//			isoer_->SetNumberOfContours(2);
//			isoer_->SetValue(0, off1);
//			isoer_->SetValue(1, off2);
//			isoer_->Update();
//			auto filler = vtkSmartPointer<vtkContourTriangulator>::New();
//			filler->SetInputConnection(isoer_->GetOutputPort());
//			filler->Update();
//			auto pprojection = vtkSmartPointer<vtkTransform>::New();
//			pprojection->Translate(0.0, 0.0, z);
//			auto transform_filter = vtkSmartPointer<vtkTransformPolyDataFilter>::New();
//			transform_filter->SetTransform(pprojection);
//			transform_filter->SetInputData(filler->GetOutput());
//			transform_filter->Update();
//			return transform_filter->GetOutput();
//		}
//
//		vtkSmartPointer<vtkPolyData> band(double off, double z1, double z2) {
//			isoer_->SetNumberOfContours(1);
//			isoer_->SetValue(0, off);
//			isoer_->Update();
//			auto extruder = vtkSmartPointer<vtkLinearExtrusionFilter>::New();
//			extruder->SetInputConnection(isoer_->GetOutputPort());
//			extruder->SetVector(0.0, 0.0, 1.0);
//			extruder->SetScaleFactor(z2 - z1);
//			extruder->Update();
//			auto pprojection = vtkSmartPointer<vtkTransform>::New();
//			pprojection->Translate(0.0, 0.0, z1);
//			auto transform_filter = vtkSmartPointer<vtkTransformPolyDataFilter>::New();
//			transform_filter->SetTransform(pprojection);
//			transform_filter->SetInputData(extruder->GetOutput());
//			transform_filter->Update();
//			return transform_filter->GetOutput();
//		}
//
//		vtkSmartPointer<vtkPolyData> tongue_and_groove(SeamParameters params) {
//			auto appender = vtkSmartPointer<vtkAppendPolyData>::New();
//			//The isocontour extraction is super redundant here fix in future!!!!
//			double z_tongue_tip = (params.gap_depth - params.tongue_depth) / 2.0;
//			double z_teeth_back = z_tongue_tip - params.gap_depth;
//			double z_teeth_front = z_teeth_back - params.trim_depth;
//			auto flat1 = flat(params.groove_outer, z_teeth_front);
//			appender->AddInputData(flat1);
//			auto band2 = band(params.groove_outer, z_teeth_front, -z_tongue_tip);
//			appender->AddInputData(band2);
//			auto flat3 = flat(params.groove_inner, params.groove_outer, -z_tongue_tip);
//			appender->AddInputData(flat3);
//			auto band4 = band(params.groove_inner, z_teeth_back, -z_tongue_tip);
//			appender->AddInputData(band4);
//			auto flat5 = flat(params.groove_inner, z_teeth_back);
//			appender->AddInputData(flat5);
//			double innest = params.groove_inner - params.gap_radial;
//			auto flat6 = flat(innest, z_tongue_tip);
//			appender->AddInputData(flat6);
//			auto band7 = band(innest, z_tongue_tip, -z_teeth_back);
//			appender->AddInputData(band7);
//			auto flat8 = flat(innest, params.groove_outer, -z_teeth_back);
//			appender->AddInputData(flat8);
//			auto band9 = band(params.groove_outer, -z_teeth_back, -z_teeth_front);
//			appender->AddInputData(band9);
//			auto flat10 = flat(params.groove_outer, -z_teeth_front);
//			appender->AddInputData(flat10);
//			appender->Update();
//			auto cleaner = vtkSmartPointer<vtkCleanPolyData>::New();
//			cleaner->PointMergingOn();
//			cleaner->SetInputData(appender->GetOutput());
//			cleaner->Update();
//			return cleaner->GetOutput();
//			//Put all that together into one mesh and merge all the duplicated points...
//		}
//		*/
	}; // end namespace seam
}; //end namespace d3d