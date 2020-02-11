#ifndef _USE_MATH_DEFINES
#define _USE_MATH_DEFINES
#endif // !1

#include <ctime>

#include "seamDesigner.h"
#include <vtkAppendPolyData.h>
#include <vtkBooleanOperationPolyDataFilter.h>
#include <vtkCenterOfMass.h>
#include <vtkCleanPolyData.h>
#include <vtkClipPolyData.h>
#include <vtkContourTriangulator.h>
#include <vtkCutter.h>
#include <vtkLinearExtrusionFilter.h>
#include <vtkPlane.h>
#include <vtkPointData.h>
#include <vtkPolyDataNormals.h>
#include <vtkTransformPolyDataFilter.h>
#include <vtkTriangleFilter.h>
#include <vtkXMLPolyDataWriter.h>

//!!!!DEBUG
void write_polydata_debug(const vtkSmartPointer<vtkPolyData> x, const std::string filename) {
	//Write
	auto writer3 = vtkSmartPointer<vtkXMLPolyDataWriter>::New();
	writer3->SetInputData(x);
	writer3->SetFileName(filename.c_str());
	writer3->Write();
}

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

		// Return the capped surface above/below the plane
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
			appender->AddInputConnection(filler->GetOutputPort());
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
		}

		vtkSmartPointer<vtkPolyData> CrossSectioner::get_planed_cross_section(const double margin) {
			auto cs = get_cross_section();
			if (!is_planar_transform_computed_) compute_planar_transform();
			auto transformer = vtkSmartPointer<vtkTransformPolyDataFilter>::New();
			transformer->SetInputData(cs);
			transformer->SetTransform(planar_transform_);
			transformer->Update();
			return transformer->GetOutput();
		}

		//
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

		void CrossSectioner::compute_planar_transform() {
			auto cs = get_cross_section();
			if (cs->GetNumberOfPoints() == 0) {
				double id_mat[16]{1,0,0,0, /**/ 0,1,0,0, /**/ 0,0,1,0, /**/ 0,0,0,1 /**/ };
				auto transform = vtkSmartPointer<vtkTransform>::New();
				transform->PostMultiply();
				transform->SetMatrix(id_mat);
				planar_transform_ = transform;
				planar_transform_->Update();
				is_planar_transform_computed_ = true;
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
				std::array<double,16> ortho;
				if (is_tongue_) {
					vtkMath::Cross(tongue_direction_.data(), planarx, planary);
					vtkMath::Normalize(planary);
					ortho = { planarx[0],planarx[1],planarx[2],0.0, planary[0],planary[1],planary[2],0.0, tongue_direction_[0], tongue_direction_[1], tongue_direction_[2],0.0,  0.0,0.0,0.0,1.0 };

				}
				else {
					vtkMath::Cross(plane_normal_.data(), planarx, planary);
					vtkMath::Normalize(planary);
					ortho = { planarx[0],planarx[1],planarx[2],0.0, planary[0],planary[1],planary[2],0.0, plane_normal_[0], plane_normal_[1], plane_normal_[2],0.0,  0.0,0.0,0.0,1.0 };
				}
				//The ppoints are the in-plane points: planar points.
				transform->Translate(-origin[0], -origin[1], -origin[2]);				
				transform->Concatenate(ortho.data());
				planar_transform_ = transform;
				planar_transform_->Update();
				is_planar_transform_computed_ = true;
			}
		};

		//
		vtkSmartPointer<vtkTransform> CrossSectioner::get_planar_transform() {
			if (!is_planar_transform_computed_) compute_planar_transform();
			return planar_transform_;
		};
//
//		/*~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~*/
//		/* Joint Tongue and Groove Geometry Designer */
//		/*~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~*/

		SeamDesigner::SeamDesigner() = default;

		SeamDesigner::SeamDesigner(d3d::CommonMeshData & meshin, SeamDesigner::Parameters & parameters)
		{
			parameters_ = parameters;
			vtk_polydata_ = makePolydata(meshin);
			SeamDesigner(vtk_polydata_, parameters);
		};

		SeamDesigner::SeamDesigner(vtkSmartPointer<vtkPolyData> meshin, SeamDesigner::Parameters & parameters) {
			parameters_ = parameters;
			vtk_polydata_ = meshin;
			if (parameters_.use_tongue) {
				cross_sectioner_ = std::make_shared<CrossSectioner>(vtk_polydata_,
					parameters_.plane_origin,
					parameters_.plane_normal,
					parameters_.tongue_direction);
			}
			else
			{
				cross_sectioner_ = std::make_shared<CrossSectioner>(vtk_polydata_,
					parameters_.plane_origin,
					parameters_.plane_normal);
			}
			// Compute the cross section.
			auto cs = cross_sectioner_->get_planed_cross_section();
			auto planar_cs = std::make_shared<planar::CurveCollection>(cs);
			// In general may need to clean up the curves?
			// Set up the distance field for computing tongue and groove geometry.
			double inmax = 1.1 * (parameters_.groove_inner - parameters_.gap_radial);
			double outmax = 1.1 * (parameters_.groove_outer);
			// compute distance field for curve offsets
			offgrid_ = planar_cs->distance_field(inmax, outmax);
			offgrid_->GetPointData()->SetActiveScalars("offset");
			isoer_->SetInputData(offgrid_);
			isoer_->ComputeScalarsOff();
		};

		vtkSmartPointer<vtkPolyData> SeamDesigner::flat(double off, double z) {
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
	
		vtkSmartPointer<vtkPolyData> SeamDesigner::flat(double off1, double off2, double z) {
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
	
		vtkSmartPointer<vtkPolyData> SeamDesigner::band(double off, double z1, double z2) {
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
	
		void SeamDesigner::compute_groove() {
			// depth offsets
			const double z_tongue_tip = (parameters_.gap_depth - parameters_.tongue_depth) / 2.0;
			const double z_teeth_back = z_tongue_tip - parameters_.gap_depth;
			const double z_teeth_front = z_teeth_back - parameters_.trim_depth;
			//
			auto appender = vtkSmartPointer<vtkAppendPolyData>::New();
			auto flat1 = flat(parameters_.groove_outer, z_teeth_front);
			appender->AddInputData(flat1);
			auto band2 = band(parameters_.groove_outer, z_teeth_front, -z_tongue_tip);
			appender->AddInputData(band2);
			auto flat3 = flat(parameters_.groove_inner, parameters_.groove_outer, -z_tongue_tip);
			appender->AddInputData(flat3);
			auto band4 = band(parameters_.groove_inner, z_teeth_back, -z_tongue_tip);
			appender->AddInputData(band4);
			auto flat5 = flat(parameters_.groove_inner, z_teeth_back);
			appender->AddInputData(flat5);
			appender->Update();
			auto cleaner = vtkSmartPointer<vtkCleanPolyData>::New();
			cleaner->PointMergingOn();
			cleaner->SetInputData(appender->GetOutput());
			cleaner->Update();
			auto triangler = vtkSmartPointer<vtkTriangleFilter>::New();
			triangler->SetInputConnection(cleaner->GetOutputPort());
			triangler->Update();
			auto normaler = vtkSmartPointer<vtkPolyDataNormals>::New();
			normaler->SetInputConnection(triangler->GetOutputPort());
			normaler->SplittingOn();
			normaler->Update();
			write_polydata_debug(normaler->GetOutput(), "testgroove.vtp");
			groove_ = normaler->GetOutput();
		}

		void SeamDesigner::compute_tongue() {
			// depth offsets
			const double z_tongue_tip = (parameters_.gap_depth - parameters_.tongue_depth) / 2.0;
			const double z_teeth_back = z_tongue_tip - parameters_.gap_depth;
			const double z_teeth_front = z_teeth_back - parameters_.trim_depth;
			//
			auto appender = vtkSmartPointer<vtkAppendPolyData>::New();
			double innest = parameters_.groove_inner - parameters_.gap_radial;
			auto flat6 = flat(innest, z_tongue_tip);
			appender->AddInputData(flat6);
			auto band7 = band(innest, z_tongue_tip, -z_teeth_back);
			appender->AddInputData(band7);
			auto flat8 = flat(innest, parameters_.groove_outer, -z_teeth_back);
			appender->AddInputData(flat8);
			auto band9 = band(parameters_.groove_outer, -z_teeth_back, -z_teeth_front);
			appender->AddInputData(band9);
			auto flat10 = flat(parameters_.groove_outer, -z_teeth_front);
			appender->AddInputData(flat10);
			appender->Update();
			auto cleaner = vtkSmartPointer<vtkCleanPolyData>::New();
			cleaner->PointMergingOn();
			cleaner->SetInputData(appender->GetOutput());
			cleaner->Update();
			auto triangler = vtkSmartPointer<vtkTriangleFilter>::New();
			triangler->SetInputConnection(cleaner->GetOutputPort());
			triangler->Update();
			auto normaler = vtkSmartPointer<vtkPolyDataNormals>::New();
			normaler->SetInputConnection(triangler->GetOutputPort());
			normaler->SplittingOn();
			normaler->Update();
			write_polydata_debug(normaler->GetOutput(), "testtongue.vtp");
			tongue_ = normaler->GetOutput();
		}

		vtkSmartPointer<vtkPolyData> SeamDesigner::get_bottom() {
			// First compute the groove in the xy-plane
			compute_groove();
			auto inv_transform = cross_sectioner_->get_planar_transform()->GetInverse();
			auto untransformer = vtkSmartPointer<vtkTransformPolyDataFilter>::New();
			untransformer->SetTransform(inv_transform);
			untransformer->SetInputData(groove_);
			untransformer->Update();
			// Join the bottom and groove
			auto unioner = vtkSmartPointer<vtkBooleanOperationPolyDataFilter>::New();
			unioner->SetOperationToUnion();
			unioner->SetInputConnection(0, untransformer->GetOutputPort());
			const double z_bot_clip = -0.5 * parameters_.tongue_depth - parameters_.gap_depth - 0.5 * parameters_.trim_depth;
			auto bottom_base = cross_sectioner_->get_clipping(z_bot_clip, false);
			auto cleaner = vtkSmartPointer<vtkCleanPolyData>::New();
			cleaner->PointMergingOn();
			cleaner->SetInputData(bottom_base);
			cleaner->Update();
			auto triangler = vtkSmartPointer<vtkTriangleFilter>::New();
			triangler->SetInputConnection(cleaner->GetOutputPort());
			triangler->Update();
			auto normaler = vtkSmartPointer<vtkPolyDataNormals>::New();
			normaler->SetInputConnection(triangler->GetOutputPort());
			normaler->SplittingOn();
			normaler->Update();
			unioner->SetInputConnection(1, normaler->GetOutputPort());
			unioner->Update();
			auto recleaner = vtkSmartPointer<vtkCleanPolyData>::New();
			recleaner->PointMergingOn();
			recleaner->SetInputConnection(unioner->GetOutputPort());
			recleaner->Update();
			return recleaner->GetOutput();
		}

		vtkSmartPointer<vtkPolyData> SeamDesigner::get_top() {
			// First compute the tongue in the xy-plane
			compute_tongue();
			auto inv_transform = cross_sectioner_->get_planar_transform()->GetInverse();
			auto untransformer = vtkSmartPointer<vtkTransformPolyDataFilter>::New();
			untransformer->SetTransform(inv_transform);
			untransformer->SetInputData(tongue_);
			untransformer->Update();
			// Join the top and tongue
			auto unioner = vtkSmartPointer<vtkBooleanOperationPolyDataFilter>::New();
			unioner->SetOperationToUnion();
			unioner->SetInputConnection(0, untransformer->GetOutputPort());
			const double z_top_clip = 0.5 * parameters_.tongue_depth + parameters_.gap_depth + 0.5 * parameters_.trim_depth;
			auto top_base = cross_sectioner_->get_clipping(z_top_clip, true);
			auto cleaner = vtkSmartPointer<vtkCleanPolyData>::New();
			cleaner->PointMergingOn();
			cleaner->SetInputData(top_base);
			cleaner->Update();
			auto triangler = vtkSmartPointer<vtkTriangleFilter>::New();
			triangler->SetInputConnection(cleaner->GetOutputPort());
			triangler->Update();
			auto normaler = vtkSmartPointer<vtkPolyDataNormals>::New();
			normaler->SetInputConnection(triangler->GetOutputPort());
			normaler->SplittingOn();
			normaler->Update();
			unioner->SetInputConnection(1, normaler->GetOutputPort());
			unioner->Update();
			auto recleaner = vtkSmartPointer<vtkCleanPolyData>::New();
			recleaner->PointMergingOn();
			recleaner->SetInputConnection(unioner->GetOutputPort());
			recleaner->Update();
			return recleaner->GetOutput();
		}

		vtkSmartPointer<vtkPolyData> SeamDesigner::get_seam() {
			auto top = get_top();
			auto bottom = get_bottom();
			auto appender = vtkSmartPointer<vtkAppendPolyData>::New();
			appender->AddInputData(top);
			appender->AddInputData(bottom);
			appender->Update();
			return appender->GetOutput();
		}
	}; // end namespace seam
}; //end namespace d3d