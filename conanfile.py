import os
from conans import ConanFile, CMake


class SoapConan(ConanFile):
    settings = "os", "compiler", "build_type", "arch"

    requires = "boost/1.70.0@conan/stable", \
                "VTK-maxi/0.1@d3d/testing", \
                "eigen/3.3.7@conan/stable"

    generators = "cmake"

    def config_options(self):
        if self.settings.os == "Windows":
            self.options["TBB"].shared = True

    def imports(self):
        self.copy("*.lib", dst="lib", keep_path=False)
        self.copy("*.dll", dst="bin", keep_path=False)
        self.copy("*.dylib*", dst="lib", keep_path=False)
        self.copy("*.so*", dst="lib", keep_path=False)
        self.copy("*.a", dst="lib", keep_path=False)