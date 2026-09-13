#ifndef LU_H
#define LU_H

int Get_Pivots_Extents(double **A, int N, double Tol, int *P, int *P_max,
                       int *Row_nz, int *Col_nz_U, int *Col_nz_L);

void Get_Stops(double **A, int N, int *P_max,
               int **Row_ids,
               int **Col_ids_U,
               int **Col_ids_L);

void LUDecompose(double **A, int N, int *P,
                   int *Row_nz, int **Row_ids,
                   int *Col_nz_U, int **Col_ids_U,
                   double **A_L, double **A_U);

void LUPSolve(double **A, int *P, double *b, int N, double *x,
              int *Col_nz_U, int **Col_ids_U,
              int *Col_nz_L, int **Col_ids_L,
              double **A_L, double **A_U);

#endif