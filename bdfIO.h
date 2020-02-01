#ifndef MESHDATALIB_BDFIO_H
#define MESHDATALIB_BDFIO_H

#include <boost/filesystem.hpp>
#include <vector>
#include "d3derr.h"
#include "meshData.h"

namespace d3d {
namespace io {
D3D_status writeBDFFromCommonMeshData(CommonMeshData& mesh,
                            const boost::filesystem::path& path);

D3D_status readBDFToCommonMeshData(const boost::filesystem::path& meshPath,
                         CommonMeshData& mesh, int& designDomainPID);
D3D_status readBDFToCommonMeshData(const boost::filesystem::path& meshPath,
                         CommonMeshData& mesh);
}  // namespace io

}  // namespace d3d

#endif