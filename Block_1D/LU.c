#include <stdlib.h>
#include <stdio.h>
#include <math.h>
#include "LU.h"

int Get_Pivots_Extents(double **A, int N, double Tol, int *P, int *P_max,
                       int *Row_nz, int *Col_nz_U, int *Col_nz_L,
                       int **Checklist_S){
    int i, j, k, imax, *iptr; 
    double maxA, *ptr, absA;
    int row_counter, col_counter;

    // Save initial values with a 1
    for (i = 0; i < N; i++)
        for (j = 0; j < N; j++)
            if (A[i][j] != 0.0)
                Checklist_S[i][j] = 1;

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

            // pivoting checklist
            iptr = Checklist_S[i];
            Checklist_S[i] = Checklist_S[imax];
            Checklist_S[imax] = iptr;
        }

        Col_nz_U[i] = 1;
        
        row_counter = 0;
        for (j = i + 1; j < N; j++) {
            if (A[j][i] != 0.0){
                row_counter++;

                A[j][i] /= A[i][i]; // one scratch for here
                if (Checklist_S[j][i] == 0)
                    Checklist_S[j][i] = 2;

                col_counter = 0;
                Col_nz_L[j]++;
                for (k = i + 1; k < N; k++){
                    A[j][k] -= A[j][i] * A[i][k]; // multiple scratches for here
                    if (A[i][k] != 0.0){
                        col_counter++;
                        if (Checklist_S[j][k] == 0)
                            Checklist_S[j][k] = 2;
                    }
                }

                Col_nz_U[i] = col_counter + 1;
            }
        }
        Row_nz[i] = row_counter;
    }

    return 1;
}

void Get_Stops(double **A, int N, int *P_max, int *P,
               int **Row_ids,
               int **Col_ids_U,
               int **Col_ids_S1,
               int **Col_ids_S2,
               int ***Col_ids_S3,
               int **Global_Row,
               int **Col_ids_L){
    int i, j, k, imax; 
    double *ptr;
    int row_counter, col_counter;
    int *P_Global = malloc(N*sizeof(int));

    for (i = 0; i < N; i++)
        P_Global[i] = i; // For global row tracking

    for (i = 0; i < N; i++) {

        imax = P_max[i];
        if (imax > i) {

            ptr = A[i];
            A[i] = A[imax];
            A[imax] = ptr;

            j = P_Global[i];
            P_Global[i] = P_Global[imax];
            P_Global[imax] = j;
        }
        // Need each rows diagonal column index, as well as the column index for the A[j][i] call
        Col_ids_U[i][0] = i;

        row_counter = 0; 
        for (j = i + 1; j < N; j++) {
            if (A[j][i] != 0.0){
                A[j][i] /= A[i][i];
                Row_ids[i][row_counter] = j;
                Col_ids_S2[i][row_counter] = i;
                Global_Row[i][row_counter] = P_Global[j];
                row_counter++;
                // Need to save off what global row, row id belongs to.
                // Need to have a separate Col_ids_U mapped to S for each row j
                // So for each i, for each j, have col mapping of k
                col_counter = 1;
                Col_ids_S3[i][row_counter-1][0] = i;
                for (k = i + 1; k < N; k++){
                    A[j][k] -= A[j][i] * A[i][k];
                    if (A[i][k] != 0.0){
                        Col_ids_S3[i][row_counter-1][col_counter] = k;
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

    free(P_Global);
}

void LUDecompose(double **A, double **A_S, int N, int *P,
                   int *Row_nz, int **Row_ids,
                   int *Col_nz_U, int **Col_ids_U,
                   int **Col_ids_S1, int **Col_ids_S2, int ***Col_ids_S3,
                   double **A_L, double **A_U) {
    int i, j, k, imax; 
    double *ptr;
    int row_id, col_id1, col_id3;
    int *A_L_counter = calloc(N, sizeof(int));

    for (i = 0; i < N; i++) {

        imax = P[i];
        if (imax > i) {

            //pivoting rows of A
            ptr = A[i];
            A[i] = A[imax];
            A[imax] = ptr;

            //pivoting rows of A_S
            ptr = A_S[i];
            A_S[i] = A_S[imax];
            A_S[imax] = ptr;

            //pivoting rows of A_L
            ptr = A_L[i];
            A_L[i] = A_L[imax];
            A_L[imax] = ptr;

            //pivoting counter
            j = A_L_counter[i];
            A_L_counter[i] = A_L_counter[imax];
            A_L_counter[imax] = j;

        }

        for (j = 0; j < Col_nz_U[i]; j++){ // This can use the same length as Col_nz_U[i], just need A_S specific column indices
            // A_U[i][j] = A[i][Col_ids_U[i][j]];
            A_U[i][j] = A_S[i][Col_ids_S1[i][j]];
        }
        
        // Here Row id is the relative row below row i, it may change from the final layout in P
        // Need to know what global row, row id is 
        for (j = 0; j < Row_nz[i]; j++) {
            row_id = Row_ids[i][j];
            col_id1 = A_L_counter[row_id];

            // A_L[row_id][col_id1] = A[row_id][i] / A_U[i][0]; // This is actually the same as the Row_nz, just need A_S specific column indices
            A_L[row_id][col_id1] = A_S[row_id][Col_ids_S3[i][j][0]] / A_U[i][0];
            A_L_counter[row_id]++;
            for (k = 1; k < Col_nz_U[i]; k++){
                // col_id2 = Col_ids_U[i][k];
                col_id3 = Col_ids_S3[i][j][k];
                // A[row_id][col_id2] -= A_L[row_id][col_id1] * A_U[i][k]; // This is the same as Col_nz_U[i], can reuse A_S specific column indices from above while starting at i+1
                A_S[row_id][col_id3] -= A_L[row_id][col_id1] * A_U[i][k];
                // if (A[row_id][col_id2] != A_S[row_id][col_id3])
                    // printf("HERE");
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
            // A[i][col_id] = 0.0;
        }
    }
    // backward substitution using upper matrix (includes diagonal)
    for (int i = N - 1; i >= 0; i--) {
        for (int k = 1; k < Col_nz_U[i]; k++){
            col_id = Col_ids_U[i][k];
            x[i] -= A_U[i][k] * x[col_id];
            // A[i][col_id] = 0.0;
        }
        x[i] /= A_U[i][0];
        // A[i][i] = 0.0;
    }
}