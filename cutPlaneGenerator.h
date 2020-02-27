#ifndef SOAP_CUT_PLANE_GENERATOR_H
#define SOAP_CUT_PLANE_GENERATOR_H

#include <Eigen/Dense>

#include <vtkCellLocator.h>
#include <vtkContourFilter.h>
#include <vtkPlane.h>
#include <vtkPolyData.h>
#include <vtkSmartPointer.h>

namespace soap {
namespace cuts {

class Plane {
	public:
		const Eigen::Vector3d origin_ {0.0,0.0,0.0};
		const Eigen::Vector3d normal_ {0.0,0.0,1.0};
		vtkSmartPointer<vtkPlane> get_vtkPlane();
		Eigen::Vector3d project_point(Eigen::Vector3d x);
		Plane() {};
		Plane(const Eigen::Vector3d origin, const Eigen::Vector3d normal);
};

class IsocurvePlanarizer {
	/* From a polydata with a scalar field attached,
	extract isocurve from the field and fit a plane.
	The plane normal is best fit orthogonal to the isocurves's surface normals,
	and the plane origin is the isocurve centroid.
	*/
private:
	vtkSmartPointer<vtkPolyData> polydata_ = vtkSmartPointer<vtkPolyData>::New();
	vtkSmartPointer<vtkContourFilter> iso_ = vtkSmartPointer<vtkContourFilter>::New();
	vtkSmartPointer<vtkCellLocator> locator_ = vtkSmartPointer<vtkCellLocator>::New();

public:
	IsocurvePlanarizer() {};
	IsocurvePlanarizer(vtkSmartPointer<vtkPolyData> polydata, const std::string scalar_array_name);
	Plane plane_estimate(vtkSmartPointer<vtkPolyData> curve);
	vtkSmartPointer<vtkPolyData> get_isocurve(const double isocurve_value);
	//std::vector<Plane> get_planes(const double isocurve_value);
};

}; //end namespace cut
};//end namsespace soap

#endif