/*----------------------------------------------------------------------*/
/*              Copyright 2018 Divergent Technologies, Inc              */
/*----------------------------------------------------------------------*/
#include <system_error>

//  D3D error code definitions

#ifndef D3DERR_CLASSDECL_H__58B47312F1664986920E84FD08E6BDF7__INCLUDED_
#define D3DERR_CLASSDECL_H__58B47312F1664986920E84FD08E6BDF7__INCLUDED_

/// General
//#define D3D_SUCCESS 0

enum class D3D_status {
    SUCCESS = 0,
    FAIL,
	NO_OPERATION,
    UNKNOWN_OPTIONS,
    MISSING_FILE_EXTENSION,
    MISSING_INPUT,
    CANNOT_READ_MESH,
    CANNOT_OPEN_FILE,
    CANNOT_WRITE_FILE,
    ERROR_INVALID_MESH_INPUT,
    ERROR_COLINEAR_TRIANGLE,
    ERROR_TETGEN_OUT_OF_MEMORY,
    ERROR_TETGEN_UNKNOWN_FAILURE,
    ERROR_TETGEN_SELF_INTERSECT,
    ERROR_TETGEN_SMALL_FEATURE,
    ERROR_TETGEN_CLOSE_FACETS,
	ERROR_IN_COMPARISON,
    FAIL_CLPS_DUP_EDGE,
    FAIL_CLPS_FLIP_NORMAL,
    FAIL_CLSP_SHAPE_DECREASE,
    FAIL_NON_MANIFOLD,
    FAIL_PATCH_HOLES,
    FAIL_SMALL_DIMENSION,
    FAIL_COTWEIGHT,
    FAIL_CLPS_BRIDGE_GEDGE,
    FAIL_SWAP_DUP_EDGE,
	FAIL_NORM_TOO_SMALL,
    INFO_NO_NON_DESIGN_SPACE,
    INFO_NO_DESIGN_SPACE_BOUNDARY,
	INFO_FOUND_CONSTRAINED_NODES,
	ERROR_NUMERICAL_TOLERANCE,
    PID_NOT_FOUND
};

namespace std {
template <>
struct is_error_code_enum<D3D_status> : public true_type {};
}  // namespace std

std::error_code make_error_code(D3D_status e);

int D3D_return(D3D_status status_code);
#endif
