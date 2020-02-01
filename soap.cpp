#include <iostream>

#include "bdfIO.h"
#include "seamDesigner.h"

int main() {
	std::cout << "Saluton Mundo!";
	d3d::CommonMeshData meshin;
	boost::filesystem::path filein = "C:\\Users\\sscott\\Pictures\\unitsphere_meshlab.bdf";
	auto err = d3d::io::readBDFToCommonMeshData(filein, meshin);
	d3d::soap::Designer::Parameters paramers;
	paramers.plane_origin = { 462.15894026125324, 119.3025104782025, -100.30062516465588 };
	paramers.plane_normal = { -0.1398971605146768, 0.9816945964030201, -0.12924590466642394 };
	paramers.tongue_direction = { 0.0, 1.0, 0.0 };
	paramers.groove_outer = 0.01;
	paramers.groove_inner = -0.01;
	paramers.gap_radial = 0.005;
	paramers.gap_depth = 0.005;
	paramers.tongue_depth = 0.05;
	paramers.trim_depth = 0.025;
	paramers.use_tongue = true;
	d3d::soap::Designer dessy = d3d::soap::Designer(meshin, paramers);
};