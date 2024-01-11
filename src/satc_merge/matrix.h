#pragma once

#include <cinttypes>
#include <cstring>
#include <algorithm>
#include <numeric>
#include <map>

//#define REFRESH_MATRIX_USE_OPENBLAS

//#define REFRESH_MATRIX_ALLOW_PRINT

#ifdef REFRESH_MATRIX_USE_OPENBLAS
#include <complex>
#define LAPACK_COMPLEX_CPP
#include <cblas.h>
#include <lapacke_config.h>
#include <lapacke.h>
#endif

namespace refresh
{
	const unsigned matrix_row_major = 0;
	const unsigned matrix_col_major = 1;
}

#include "matrix_1d.h"
#include "matrix_dense.h"
#include "matrix_sparse.h"
#include "matrix_sparse_compact.h"
