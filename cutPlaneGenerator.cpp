#include "cutPlaneGenerator.h"

#include <vtkCellData.h>
#include <vtkGenericCell.h>
#include <vtkPointData.h>
#include <vtkPolyDataNormals.h>

namespace soap {
namespace cuts {

	template <typename T>
	class Averager {
	/*For computing average vectors and matrices.*/
	private:
		T value;
		int num_samples{ 0 };
	public:
		Averager() {};
		Averager(T x) { value = x; }
		void include_sample(T x) {
			value = num_samples * value + x;
			num_samples++;
			value /= num_samples;
		}
		T get_value() {
			return value;
		};
	};

	vtkSmartPointer<vtkPlane> Plane::get_vtkPlane() {
		auto a_vtk_plane = vtkSmartPointer<vtkPlane>::New();
		a_vtk_plane->SetOrigin(origin_.data());
		a_vtk_plane->SetNormal(normal_.data());
		return a_vtk_plane;
	}

	Plane::Plane(const Eigen::Vector3d origin, const Eigen::Vector3d normal) :
		origin_(origin),
		normal_(normal)
	{};

	Eigen::Vector3d Plane::project_point(Eigen::Vector3d x) {
		Eigen::Vector3d pp(x - origin_);
		return origin_ + pp - pp.dot(normal_) * normal_;
	}

	IsocurvePlanarizer::IsocurvePlanarizer(vtkSmartPointer<vtkPolyData> polydata, const std::string scalar_array_name)
	{
		auto normaler = vtkSmartPointer<vtkPolyDataNormals>::New();
		normaler->SetInputData(polydata);
		normaler->ComputeCellNormalsOn();
		normaler->ComputePointNormalsOff();
		normaler->Update();
		polydata_ = normaler->GetOutput();
		polydata_->GetPointData()->SetActiveScalars(scalar_array_name.c_str());
		iso_->SetInputData(polydata_);
		locator_->SetDataSet(polydata_);
		locator_->BuildLocator();
	}
	
	Plane IsocurvePlanarizer::plane_estimate(vtkSmartPointer<vtkPolyData> curve) {
		/*
		Estimate the plane whose origin is the curve centroid and normal most orthogonal to the surface normals along curve.
		*/
		Averager<Eigen::Vector3d> centroid_avger;
		Eigen::VectorXd r6origin(6);
		r6origin << 0, 0, 0, 0, 0, 0;
		Averager<Eigen::VectorXd> nrm_quad_form_avger(r6origin);
		//auto lines = curve->GetLines();
		auto _gencell = vtkSmartPointer<vtkGenericCell>::New();
		for (int foo = 0; foo < curve->GetNumberOfLines(); foo++) {
			auto pt_list = vtkSmartPointer<vtkIdList>::New();
			curve->GetCellPoints(foo, pt_list);
			int pid0 = pt_list->GetId(0);
			int pid1 = pt_list->GetId(1);
			double pt_c0[3], pt_c1[3];
			curve->GetPoint(pid0, pt_c0);
			curve->GetPoint(pid1, pt_c1);
			Eigen::Vector3d pt0(pt_c0[0],pt_c0[1], pt_c0[2]);
			Eigen::Vector3d pt1(pt_c1[0], pt_c1[1], pt_c1[2]);
			Eigen::Vector3d midp = 0.5*(pt0+pt1);
			centroid_avger.include_sample(midp);
			//Get the normal of the closeset cell to this point.
			double _c[3], _d;
			vtkIdType cellid;
			int _b;
			locator_->FindClosestPoint(midp.data(), _c, _gencell, cellid, _b, _d);
			double nrm[3];
			polydata_->GetCellData()->GetNormals()->GetTuple(cellid, nrm);
			//Eigen::Vector3d nrm(nrm_cid);
			Eigen::VectorXd nrmsq(6);
			nrmsq << (nrm[0] * nrm[0]), (nrm[0] * nrm[1]), (nrm[0] * nrm[2]), (nrm[1] * nrm[1]), (nrm[1] * nrm[2]), (nrm[2] * nrm[2]);
			nrm_quad_form_avger.include_sample(nrmsq);
		}
		// The normal direction of the plane is the eigenvector for the least eigenvalue
		// of the quadratic form sum_{n normal vector} |n><n|
		Eigen::Vector3d origin = centroid_avger.get_value();
		Eigen::Matrix3d nrm_quad_form;
		auto nqf = nrm_quad_form_avger.get_value();
		nrm_quad_form << nqf[0], nqf[1], nqf[2],
						 nqf[1], nqf[3], nqf[4],
						 nqf[2], nqf[4], nqf[5];
		Eigen::SelfAdjointEigenSolver<Eigen::Matrix3d> eiger(nrm_quad_form);
		auto least_vec = eiger.eigenvectors().col(0);
		std::cout << origin.transpose() << std::endl;
		std::cout << least_vec.transpose() << std::endl;
		Plane aplane(origin, least_vec);
		return aplane;
	}
	vtkSmartPointer<vtkPolyData> IsocurvePlanarizer::get_isocurve(const double isocurve_value)
	{
		iso_->SetValue(0, isocurve_value);
		iso_->Update();
		return iso_->GetOutput();
	};

	//std::vector<Plane> IsocurvePlanarizer::get_planes(const double isocurve_value)
	//{
	//	iso_->SetValue(0, isocurve_value);
	//	iso_->Update();
	//	auto isocurve = iso_->GetOutput();
	//}
	

}; //end namespace cuts
}; //end namespace soap