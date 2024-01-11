#include "ob_utils.h"
#include <cinttypes>
#include <algorithm>
#include <iostream>

using namespace std;

bool ob_svd(refresh::matrix_dense<double, refresh::matrix_col_major>& svdmat, refresh::matrix_1d<double>& u0, refresh::matrix_1d<double>& v0)
{
    int64_t n_rows = svdmat.rows();
    int64_t n_cols = svdmat.cols();
    int64_t n_min = std::min(n_rows, n_cols);

    refresh::matrix_dense<double, refresh::matrix_col_major> U(n_rows, n_min);
    refresh::matrix_dense<double, refresh::matrix_col_major> Vt(n_min, n_cols);
    refresh::matrix_1d<double> S(n_min);

    lapack_int res_sdd = LAPACKE_dgesdd(
        LAPACK_COL_MAJOR, 'S',
        n_rows, n_cols,
        svdmat.data(), n_rows,
        S.data(),
        U.data(), n_rows,
        Vt.data(),
        n_min);

    u0 = U.get_col(0);
    v0 = Vt.get_row(0);

    return true;
}
