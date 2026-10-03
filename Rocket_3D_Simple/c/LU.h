#ifndef LU_H
#define LU_H

int Get_Pivots_Extents(double **A, int N, double Tol, int *P, int *P_max,
                       int *Row_nz, int *Col_nz_U, int *Col_nz_L,
                       int **Checklist_S);

void Get_Stops(double **A, int N, int *P_max, int *P,
               int **Row_ids,
               int **Col_ids_U,
               int **Col_ids_S1,
               int **Col_ids_S2,
               int ***Col_ids_S3,
               int **Global_Row,
               int **Col_ids_L);

void LUDecompose(double **A_S, int N, int *P,
                 int *Row_nz, int **Row_ids, int *Col_nz_U,
                 int **Col_ids_S1, int ***Col_ids_S3);

void LUPSolve(double **A_S, int *P, double *b, int N, double *x,
              int *Col_nz_U, int **Col_ids_U,
              int *Col_nz_L, int **Col_ids_L,
              int **Col_ids_S1, int **Col_ids_S2);

#endif