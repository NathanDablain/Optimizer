#ifndef LU_H
#define LU_H

int Get_Pivots_Extents(double **A, int N, double Tol, int *P, int *P_max, int *Row_nz, int *Col_nz, int *Col_nz_fw);

void Get_Stops(double **A, int N, int *P_max, int *Row_nz, int *Col_nz, int *Col_nz_fw, int **Stop_rows, int **Stop_cols, int **Stop_cols_fw);

void LUPDecompose2(double **A, int N, int *P, int *Row_nz, int *Col_nz, int **Stop_rows, int **Stop_cols);

void LUPSolve(double **A, int *P, double *b, int N, double *x, int *Col_nz, int **Stop_cols, int *Col_nz_fw, int **Stop_cols_fw);

#endif