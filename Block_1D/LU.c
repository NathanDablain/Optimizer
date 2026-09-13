#include <stdlib.h>
#include <stdio.h>
#include <math.h>
#include "LU.h"

int Get_Pivots_Extents(double **A, int N, double Tol, int *P, int *P_max,
                       int *Row_nz, int *Col_nz_U, int *Col_nz_L){
    int i, j, k, imax; 
    double maxA, *ptr, absA;
    int row_counter, col_counter;

    for (i = 0; i < N; i++)
        P[i] = i; //Unit permutation matrix, P[N] initialized with N

    for (i = 0; i < N; i++) {
        maxA = 0.0;
        imax = i;

        for (k = i; k < N; k++)
            if ((absA = fabs(A[k][i])) > maxA) { 
                maxA = absA;
                imax = k;
            }

        if (maxA < Tol) return 0; //failure, matrix is degenerate
        
        P_max[i] = imax;
        if (imax != i) {
            //pivoting P
            j = P[i];
            P[i] = P[imax];
            P[imax] = j;

            //pivoting rows of A
            ptr = A[i];
            A[i] = A[imax];
            A[imax] = ptr;

            //pivoting A_L
            j = Col_nz_L[i];
            Col_nz_L[i] = Col_nz_L[imax];
            Col_nz_L[imax] = j;
        }

        row_counter = 0;
        for (j = i + 1; j < N; j++) {
            if (A[j][i] != 0.0){
                row_counter++;

                A[j][i] /= A[i][i];
                col_counter = 0;
                Col_nz_L[j]++;
                for (k = i + 1; k < N; k++){
                    A[j][k] -= A[j][i] * A[i][k];
                    if (A[i][k] != 0.0){
                        col_counter++;
                    }
                }
                Col_nz_U[i] = col_counter;
            }
        }
        Row_nz[i] = row_counter;
    }

    return 1;  //decomposition done 
}

void Get_Stops(double **A, int N, int *P_max,
               int **Row_ids,
               int **Col_ids_U,
               int **Col_ids_L){
    int i, j, k, imax; 
    double *ptr;
    int row_counter, col_counter;

    for (i = 0; i < N-1; i++) {

        imax = P_max[i];
        if (imax > i) {

            // Need to be able to swap rows of A easily
            //pivoting rows of A
            ptr = A[i];
            A[i] = A[imax];
            A[imax] = ptr;

        }

        // Need to be able to iterate down a column and stop on non zeros
        // Alternatively, save off rows to jump to when eliminating on a certain row
        // This would look like an int double pointer
        row_counter = 0; 
        for (j = i + 1; j < N; j++) {
            if (A[j][i] != 0.0){
                A[j][i] /= A[i][i];
                Row_ids[i][row_counter] = j;
                row_counter++;
                // Need to be able to iterate along a row and stop on the relevant columns
                // Row i will be in its final form , row j can change and may have new entries added or taken away
                // Although all that really matters is what columns in row i remain which can be pulled from final values
                col_counter = 0;
                for (k = i + 1; k < N; k++){
                    A[j][k] -= A[j][i] * A[i][k];
                    if (A[i][k] != 0.0){
                        Col_ids_U[i][col_counter] = k;
                        col_counter++;
                    }
                }
            }
        }
    }

    for (i = 0; i < N; i++){
        col_counter = 0;
        for (j = 0; j < i; j++){
            if (A[i][j] != 0.0){
                Col_ids_L[i][col_counter] = j;
                col_counter++;
            }
        }
    }

}

void LUDecompose(double **A, int N, int *P,
                   int *Row_nz, int **Row_ids,
                   int *Col_nz_U, int **Col_ids_U,
                   double **A_L, double **A_U) {
    int i, j, k, imax; 
    double *ptr;
    int row_id, col_id1, col_id2;
    int *A_L_counter = calloc(N, sizeof(int));

    for (i = 0; i < N; i++) {

        imax = P[i];
        if (imax > i) {

            //pivoting rows of A
            ptr = A[i];
            A[i] = A[imax];
            A[imax] = ptr;

            //pivoting rows of A_L
            ptr = A_L[i];
            A_L[i] = A_L[imax];
            A_L[imax] = ptr;

            //pivoting counter
            j = A_L_counter[i];
            A_L_counter[i] = A_L_counter[imax];
            A_L_counter[imax] = j;
        }

        A_U[i][0] = A[i][i];
        for (j = 1; j < Col_nz_U[i]+1; j++)
            A_U[i][j] = A[i][Col_ids_U[i][j-1]];

        // The only difficulty is how to access column i?
        for (j = 0; j < Row_nz[i]; j++) {
            row_id = Row_ids[i][j];
            col_id1 = A_L_counter[row_id];

            A_L[row_id][col_id1] = A[row_id][i] / A_U[i][0];
            A_L_counter[row_id]++;
            for (k = 0; k < Col_nz_U[i]; k++){
                col_id2 = Col_ids_U[i][k];
                A[row_id][col_id2] -= A_L[row_id][col_id1] * A_U[i][k+1];

            }
        }
    }

    free(A_L_counter);
}

/* INPUT: A,P filled in LUPDecompose; b - rhs vector; N - dimension
 * OUTPUT: x - solution vector of A*x=b
 */
void LUPSolve(double **A, int *P, double *b, int N, double *x,
              int *Col_nz_U, int **Col_ids_U,
              int *Col_nz_L, int **Col_ids_L,
              double **A_L, double **A_U) {
    int col_id;
    // forward substitution using lower matrix
    for (int i = 0; i < N; i++) {
        x[i] = b[P[i]];

        for (int k = 0; k < Col_nz_L[i]; k++){
            col_id = Col_ids_L[i][k];
            x[i] -= A_L[i][k] * x[col_id];
            A[i][col_id] = 0.0;
        }
    }
    // backward substitution using upper matrix (includes diagonal)
    for (int i = N - 1; i >= 0; i--) {
        for (int k = 0; k < Col_nz_U[i]; k++){
            col_id = Col_ids_U[i][k];
            x[i] -= A_U[i][k+1] * x[col_id];
            A[i][col_id] = 0.0;
        }
        x[i] /= A_U[i][0];
        A[i][i] = 0.0;
    }
}
