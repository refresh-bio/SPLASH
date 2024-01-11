#pragma once

#include <complex>
#include <algorithm>

using namespace std;

#define LAPACK_COMPLEX_CPP

#define OPENBLAS_OS_WINNT

#include <openblas_config.h>
#include <cblas.h>

#define USE_OPEN_BLAS

#include <lapacke_config.h>
#include <lapacke.h>

#include "matrix.h"

bool ob_svd(refresh::matrix_dense<double, refresh::matrix_col_major>& svdmat, refresh::matrix_1d<double>& u0, refresh::matrix_1d<double>& v0);

