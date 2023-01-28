/***********************************************************************
Numerical Recipes routine: gaussj.c for Gauss-Jordan Elimination
Version 3.0, May 17, 2011.
************************************************************************/
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include "nrutil.h"
//#define SWAP(a,b) {double temp=(a);(a)=(b);(b)=temp;}

void print_vector_double_isp(double **a, int n, const char *msg, int col);
void print_matrix_double(double **a, int n, const char *msg);

void gaussj(double **a, int n, double **b, int m)
{
    for (int i = 1; i <= n; i++)
    {
        // Get maximum index.
        int i_max = i;
        double maximum = a[i][i];
        for (int row = i + 1; row <= n; row++)
        {
            if (a[row][i] > maximum)
            {
                i_max = row;
                maximum = a[row][i];
            }
        }

        // Swap row i with row i_max.
        for (int col = 1; col <= n; col++)
        {
            double temp = a[i][col];
            a[i][col] = a[i_max][col];
            a[i_max][col] = temp;
        }
        for (int col = 1; col <= m; col++)
        {
            double temp = b[i][col];
            b[i][col] = b[i_max][col];
            b[i_max][col] = temp;
        }

        // Divide by diagonal entry.
        double diag = a[i][i];
        for (int col = 1; col <= n; col++) a[i][col] /= diag;
        for (int col = 1; col <= m; col++) b[i][col] /= diag;

        // Perform row linear combinations.
        for (int row = 1; row <= n; row++)
        {
            if (row != i)
            {
                double val = a[row][i];
                for (int col = 1; col <= n; col++) a[row][col] -= val * a[i][col];
                for (int col = 1; col <= m; col++) b[row][col] -= val * b[i][col];
            }
        }

//        printf("i = %d\n\n", i);
//        print_matrix_double(a, n, "a");
//        print_vector_double_isp(b, n, "b", 1);
    }
}

//void gaussj(double **a, int n, double **b, int m)
//{
//    for (int i = 1; i <= n; i++)
//    {
//        // Get argmax row.
//        int i_max = i;
//        double max = a[i][i];
//        for (int j = i + 1; j <= n; j++)
//        {
//            if (a[i][j] > max)
//            {
//                i_max = j;
//                max = a[i][j];
//            }
//        }
//        printf("imax = %d\n\n", i_max);
//
//        // Swap rows, so that row i is row i_max.
//        for (int j = 1; j <= n; j++)
//        {
//            double temp = a[i][j];
//            a[i][j] = a[i_max][j];
//            a[i_max][j] = temp;
//        }
//        for (int j = 1; j <= m; j++)
//        {
//            double temp = b[i][j];
//            b[i][j] = b[i_max][j];
//            b[i_max][j] = temp;
//        }
//
//        // Divide row i by diagonal entry.
//        for (int j = 1; j <= n; j++) if (j != i) a[i][j] /= a[i][i];
//        for (int j = 1; j <= m; j++) b[i][j] /= a[i][i];
//        a[i][i] = 1.0;
//
//        // Perform row linear combinations.
//        for (int row = 1; row <= n; row++) if (row != i)
//            {
//                double val = a[row][i];
//                for (int col = 1; col <= n; col++)
//                    if (col != i)
//                        a[row][col] -= val * a[i][col];
//                for (int col = 1; col <= m; col++)
//                    b[row][col] -= val * b[i][col];
//
//                a[row][i] = 0.0;
//            }
//
//        printf("i = %d\n\n", i);
//        print_matrix_double(a, n, "a");
//        print_vector_double_isp(b, n, "b", 1);
//    }
//}

//void gaussj(double **a, int n, double **b, int m)
//{
//    int *indxc, *indxr, *ipiv;
//    int i, icol, irow, j, k, l, ll;
//    double big, dum, pivinv;
//    indxc = ivector(1, n);
//    indxr = ivector(1, n);
//    ipiv = ivector(1, n);
//    for (j = 1; j <= n; j++) ipiv[j] = 0;
//    for (i = 1; i <= n; i++) {
//        big = 0.0;
//        for (j = 1; j <= n; j++) if (ipiv[j] != 1) for (k = 1; k <= n; k++) if (ipiv[k] == 0) {
//                        if (fabs(a[j][k]) >= big) {
//                            big = fabs(a[j][k]);
//                            irow = j;
//                            icol = k;
//                        }
//                    }
//                    else if (ipiv[k] > 1) nrerror("GAUSSJ: Singular Matrix-1");
//        ++(ipiv[icol]);
//        if (irow != icol) {
//            for (l = 1; l <= n; l++) SWAP(a[irow][l], a[icol][l])
//            for (l = 1; l <= m; l++) SWAP(b[irow][l], b[icol][l])
//        }
//        indxr[i] = irow;
//        indxc[i] = icol;
//        if (a[icol][icol] == 0.0) nrerror("GAUSSJ: Singular Matrix-2");
//        pivinv = 1.0 / a[icol][icol];
//        a[icol][icol] = 1.0;
////        for (l = 1; l <= n; l++) a[icol][l] *= pivinv;
//        for (l = 1; l <= n; l++) if (l != icol) a[icol][l] *= pivinv;
//        for (l = 1; l <= m; l++) b[icol][l] *= pivinv;
//        for (ll = 1; ll <= n; ll++) if (ll != icol) {
//                dum = a[ll][icol];
//                a[ll][icol] = 0.0;
////                for (l = 1; l <= n; l++) a[ll][l] -= a[icol][l] * dum;
////                for (l = 1; l <= m; l++) b[ll][l] -= b[icol][l] * dum;
//                for (l = 1; l <= n; l++) if (l != icol) a[ll][l] -= a[icol][l] * dum;
//                for (l = 1; l <= m; l++) if (l != icol) b[ll][l] -= b[icol][l] * dum;
//            }
//    }
//    for (l = n; l >= 1; l--) if (indxr[l] != indxc[l]) for (k = 1; k <= n; k++) SWAP(a[k][indxr[l]], a[k][indxc[l]]);
//    free_ivector(ipiv, 1, n);
//    free_ivector(indxr, 1, n);
//    free_ivector(indxc, 1, n);
//}
//#undef SWAP