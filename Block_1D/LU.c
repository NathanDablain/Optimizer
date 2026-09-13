#include <math.h>
#include "LU.h"

int Get_Pivots_Extents(double **A, int N, double Tol, int *P, int *P_max, int *Row_nz, int *Col_nz, int *Col_nz_fw){
    int i, j, k, imax; 
    double maxA, *ptr, absA;
    int row_counter, col_counter;

    for (i = 0; i <= N; i++)
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

        }

        row_counter = 0;
        for (j = i + 1; j < N; j++) {
            if (A[j][i] != 0.0){
                row_counter++;

                A[j][i] /= A[i][i];
                col_counter = 0;
                for (k = i + 1; k < N; k++){
                    A[j][k] -= A[j][i] * A[i][k];
                    if (A[i][k] != 0.0){
                        col_counter++;
                    }
                }
                Col_nz[i] = col_counter;
            }
        }
        Row_nz[i] = row_counter;
    }

    for (i = 0; i < N; i++){
        col_counter = 0;
        for (j = 0; j < i; j++){
            if (A[i][j] != 0.0)
                col_counter++;
        }
        Col_nz_fw[i] = col_counter;
    }

    return 1;  //decomposition done 
}

void Get_Stops(double **A, int N, int *P_max, int *Row_nz, int *Col_nz, int *Col_nz_fw, int **Stop_rows, int **Stop_cols, int **Stop_cols_fw){
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
                Stop_rows[i][row_counter] = j;
                row_counter++;
                // Need to be able to iterate along a row and stop on the relevant columns
                // Row i will be in its final form , row j can change and may have new entries added or taken away
                // Although all that really matters is what columns in row i remain which can be pulled from final values
                col_counter = 0;
                for (k = i + 1; k < N; k++){
                    A[j][k] -= A[j][i] * A[i][k];
                    if (A[i][k] != 0.0){
                        Stop_cols[i][col_counter] = k;
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
                Stop_cols_fw[i][col_counter] = j;
                col_counter++;
            }
        }
    }

}

void LUPDecompose2(double **A, int N, int *P, int *Row_nz, int *Col_nz, int **Stop_rows, int **Stop_cols){
    int i, j, k, imax; 
    double *ptr;
    int row_id, col_id;

    for (i = 0; i < N-1; i++) {

        imax = P[i];
        if (imax > i) {

            //pivoting rows of A
            ptr = A[i];
            A[i] = A[imax];
            A[imax] = ptr;

        }

        for (j = 0; j < Row_nz[i]; j++) {
            row_id = Stop_rows[i][j];
            A[row_id][i] /= A[i][i];

            for (k = 0; k < Col_nz[i]; k++){
                col_id = Stop_cols[i][k];
                A[row_id][col_id] -= A[row_id][i] * A[i][col_id];
            }
        }
    }
}

/* INPUT: A,P filled in LUPDecompose; b - rhs vector; N - dimension
 * OUTPUT: x - solution vector of A*x=b
 */
void LUPSolve(double **A, int *P, double *b, int N, double *x, int *Col_nz, int **Stop_cols, int *Col_nz_fw, int **Stop_cols_fw) {
    int col_id;

    for (int i = 0; i < N; i++) {
        x[i] = b[P[i]];

        for (int k = 0; k < Col_nz_fw[i]; k++){
            col_id = Stop_cols_fw[i][k];
            x[i] -= A[i][col_id] * x[col_id];
            A[i][col_id] = 0.0;
        }
    }

    for (int i = N - 1; i >= 0; i--) {
        for (int k = 0; k < Col_nz[i]; k++){
            col_id = Stop_cols[i][k];
            x[i] -= A[i][col_id] * x[col_id];
            A[i][col_id] = 0.0;
        }
        x[i] /= A[i][i];
        A[i][i] = 0.0;
    }
}
