#include <math.h>
#include <stdlib.h>
#include <stdio.h>
#include <stdbool.h>
#include <string.h>
#include "SQP.h"
#include "Block_1D.h"
#include "LU.h"

double **A;
double b[SYSTEM_SIZE] = {0};
double x[SYSTEM_SIZE] = {0};
double z[N_DECISION_VARIABLES] = {0};
double lambda[N_CONSTRAINTS] = {0};
double grad[N_DECISION_VARIABLES] = {0};
double eq[N_CONSTRAINTS] = {0};
double t[N_KNOTS] = {0};
double h[N_KNOTS-1] = {0};

// Speedup params
int P[SYSTEM_SIZE];
int P_max[SYSTEM_SIZE];

int Row_nz[SYSTEM_SIZE] = {0};
int Col_nz_U[SYSTEM_SIZE] = {0};
int Col_nz_L[SYSTEM_SIZE] = {0};
int Col_ic_S[SYSTEM_SIZE] = {0};
int Col_nz_S[SYSTEM_SIZE] = {0};
int Sparse_A_Size[SYSTEM_SIZE] = {0};
int Col_nz_J[N_DECISION_VARIABLES] = {0};

int **Row_ids;
// These are the column indices for each row above the diagonal (elimination, backwards substitution)
int **Col_ids_U;
// These are the column indices for each row below the diagonal (forwards substitution)
int **Col_ids_L;
// These are the column indices corresponding to the row ids
int **Col_ids_S1;
// These are the column indices corresponding to Col_ids_U
int **Col_ids_S2;
int ***Col_ids_S3;
int **Global_Row;
int **Zero_cols;

// Sparse matrices
double **A_S;
int **Checklist_S;
int **A_S_ids;
double **J_S;
int **J_S_ids;
int **Lambda_ids;

void Initial_Problem_Data(){
    int i;
    double tf, dt;
    
    // Fill t and h
    tf = Get_tf();
    dt = tf / (double)(N_KNOTS - 1);
    for (i = 1; i < N_KNOTS; i++){
        t[i] = t[i-1] + dt;
        h[i-1] = dt; 
    }

    // Set initial conditions
    Load_ic(z);

    // Set slack variables
    Load_slack(z);

}

void Initial_Solve(){
    int i, j;

    // Allocate A
    A = malloc(SYSTEM_SIZE*sizeof(double*));
    for (i = 0; i < SYSTEM_SIZE; i++)
        A[i] = calloc(SYSTEM_SIZE, sizeof(double));

    // Allocate Sparse A checklist
    Checklist_S = malloc(SYSTEM_SIZE*sizeof(int*));
    for (i = 0; i < SYSTEM_SIZE; i++)
        Checklist_S[i] = calloc(SYSTEM_SIZE, sizeof(int));

    // Find Pivots and Extents
    Load_Problem_Data(A, b, z, h, false);

    Get_Pivots_Extents(A, SYSTEM_SIZE, 1.0e-20, P, P_max, Row_nz, Col_nz_U, Col_nz_L, Checklist_S);
    Parse_Checklist_S();

    // Allocate sparse A
    A_S = malloc(SYSTEM_SIZE*sizeof(double*));
    for (i = 0; i < SYSTEM_SIZE; i++)
        A_S[i] = calloc(Sparse_A_Size[i], sizeof(double));

    for (i = 0; i < SYSTEM_SIZE; i++)
        memset(A[i], 0, SYSTEM_SIZE*sizeof(double));
    
    // Allocate stop indices
    Row_ids = malloc(SYSTEM_SIZE*sizeof(int*));
    Col_ids_S3 = malloc(SYSTEM_SIZE*sizeof(int**));
    Global_Row = malloc(SYSTEM_SIZE*sizeof(int*));
    for (i = 0; i < SYSTEM_SIZE; i++){
        Col_ids_S3[i] = malloc(Row_nz[i]*sizeof(int*));
        for (j = 0; j < Row_nz[i]; j++){
            Col_ids_S3[i][j] = malloc(Col_nz_U[i]*sizeof(int));
            for (int k = 0; k < Col_nz_U[i]; k++)
                Col_ids_S3[i][j][k] = -1;
        }
        Row_ids[i] = malloc(Row_nz[i]*sizeof(int));
        Global_Row[i] = malloc(Row_nz[i]*sizeof(int));
    }

    Col_ids_U = malloc(SYSTEM_SIZE*sizeof(int*));
    Col_ids_S1 = malloc(SYSTEM_SIZE*sizeof(int*));
    for (i = 0; i < SYSTEM_SIZE; i++){
        Col_ids_U[i] = malloc(Col_nz_U[i]*sizeof(int));
        Col_ids_S1[i] = malloc(Col_nz_U[i]*sizeof(int));
    }

    Col_ids_L = malloc(SYSTEM_SIZE*sizeof(int*));
    Col_ids_S2 = malloc(SYSTEM_SIZE*sizeof(int*));
    for (i = 0; i < SYSTEM_SIZE; i++){
        Col_ids_L[i] = malloc(Col_nz_L[i]*sizeof(int));
        Col_ids_S2[i] = malloc(Col_nz_L[i]*sizeof(int));
    }

    // Find stop indices
    Load_Problem_Data(A, b, z, h, false);
    Get_Stops(A, SYSTEM_SIZE, P_max, P, Row_ids, Col_ids_U, Col_ids_S1, Col_ids_S2, Col_ids_S3, Global_Row, Col_ids_L);
    Convert_id_U_S();

    Load_Problem_Data(A_S, b, z, h, true);

    Fill_A_S_Zeros();

    for (i = 0; i < SYSTEM_SIZE; i++){
        memset(A[i], 0, SYSTEM_SIZE*sizeof(double));
    }

}

void Cleanup(){
    int i, j;
    for (i = 0; i < SYSTEM_SIZE; i++){
        free(A[i]);
        free(Row_ids[i]);
        free(Col_ids_U[i]);
        free(Col_ids_L[i]);

        free(A_S[i]);
        free(Checklist_S[i]);
        free(Col_ids_S1[i]);
        free(Col_ids_S2[i]);
        for (j = 0; j < Row_nz[i]; j++)
            free(Col_ids_S3[i][j]);
        free(Col_ids_S3[i]);
        free(Zero_cols[i]);
        if (i < N_DECISION_VARIABLES){
            free(J_S_ids[i]);
            free(Lambda_ids[i]);
        }
    }
    free(A);
    free(Row_ids);
    free(Col_ids_U);
    free(Col_ids_L);

    free(A_S);
    free(Checklist_S);
    free(Col_ids_S1);
    free(Col_ids_S2);
    free(Col_ids_S3);
    free(Zero_cols);
    free(J_S_ids);
    free(Lambda_ids);
}

sqp_results Run_SQP(){
    sqp_results results;
    double delta_z[N_DECISION_VARIABLES];
    double lambda_new[N_CONSTRAINTS];
    double delta_lambda[N_CONSTRAINTS];
    const int max_iterations = 10;
    int i, iterations;
    double pf, df;
    const double pf_tolerance = 1.0e-6;
    const double df_tolerance = 1.0e-3;
    bool converged = false;

    // Load A and b matrices
    Load_Problem_Data(A_S, b, z, h, true);
    Fill_A_S_Zeros();
    for (iterations = 0; iterations < max_iterations; iterations++){

        // Setup A in form A=(L-E)+U by partial pivoting and eliminating
        LUDecompose(A_S, SYSTEM_SIZE, P_max,
                    Row_nz, Row_ids, Col_nz_U,
                    Col_ids_S1, Col_ids_S3);

        // Perform forward and backward substitution to solve for x
        LUPSolve(A_S, P, b, SYSTEM_SIZE, x,
                Col_nz_U, Col_ids_U,
                Col_nz_L, Col_ids_L,
                Col_ids_S1, Col_ids_S2);

        // A_S has been permuted, reset it so our columns line up
        Reset_A_S();

        // Extract change in decision variables and new lagrange multipliers from x
        memcpy(delta_z, x, N_DECISION_VARIABLES*sizeof(double));
        memcpy(lambda_new, &x[N_DECISION_VARIABLES], N_CONSTRAINTS*sizeof(double));
        for (i = 0; i < N_CONSTRAINTS; i++)
            delta_lambda[i] = lambda_new[i] - lambda[i];

        // Check slack variables and perform ternary search to determine what value of alpha to use
        double alpha = Select_Alpha(z, delta_z);

        // Update decision variables and lagrange multipliers with selected alpha
        for (i = 0; i < N_DECISION_VARIABLES; i++)
            z[i] += alpha*delta_z[i];
        
        for (i = 0; i < N_CONSTRAINTS; i++)
            lambda[i] += alpha*delta_lambda[i];
        
        // Load A and b matrices
        Load_Problem_Data(A_S, b, z, h, true);
        Fill_A_S_Zeros();
        
        // See if we have satisfied tolerances and can exit early
        pf = Primary_Feasability_Check();
        df = Dual_Feasability_Check();
        if (pf < pf_tolerance && df < df_tolerance){
            converged = true;
            break;
        }
    }
    results.converged = converged;
    results.iterations = iterations;
    results.pf = pf;
    results.df = df;

    return results;
}

double Primary_Feasability_Check(){
    // Checks if we are satisfying our constraints, is the result physically feasible?
    int i;
    double max_eq;

    max_eq = 0.0;
    for (i = 0; i < N_CONSTRAINTS; i++){
        eq[i] = -b[N_DECISION_VARIABLES+i];
        if (eq[i] > max_eq)
            max_eq = eq[i];
    }

    return max_eq;
}

double Dual_Feasability_Check(){
    // Checks if we are stationary
    // Find the maximum absolute value of the gradient of the lagragian 
    int i, j;
    double max_l_grad;
    double lagrange_grad[N_DECISION_VARIABLES];
    // lagrange_grad = (grad - jacobian'*lambda);

    max_l_grad = 0.0;
    for (i = 0; i < N_DECISION_VARIABLES; i++){
        grad[i] = -b[i];
        lagrange_grad[i] = grad[i];
        for (j = 0; j < Col_nz_J[i]; j++)
            lagrange_grad[i] -= A_S[i][J_S_ids[i][j]] * lambda[Lambda_ids[i][j]];

        if (fabs(lagrange_grad[i]) > max_l_grad)
            max_l_grad = fabs(lagrange_grad[i]);
    }

    return max_l_grad;
}

void Reset_A_S(){
    double *ptr;
    int i, j;
    int P_temp[SYSTEM_SIZE];
    memcpy(P_temp, P, SYSTEM_SIZE*sizeof(int));

    for (i = 0; i < SYSTEM_SIZE; i++){
        while(P_temp[i] != i){
            j = P_temp[i];
            P_temp[i] = P_temp[j];
            P_temp[j] = j;

            ptr = A_S[i];
            A_S[i] = A_S[j];
            A_S[j] = ptr;
        }
    }
}

int Search_For_Sparse_Column(int row, int col){
    int i;
    int result = -1;
    for (i = 0; i < Sparse_A_Size[row]; i++)
        if (A_S_ids[row][i] == col)
            return i;

    return result;
}

void Get_IC_Cols(){
    int i, offset_z, offset_c;

    offset_z = 0;
    offset_c = 0;
    Load_First_Knot_Columns(offset_z, offset_c);
    offset_c += FIRST_KNOT_CONSTRAINTS;

    for (i = 1; i < N_KNOTS - 1; i++){
        offset_z = i*KNOT_SIZE;
        Load_Middle_Knot_Columns(i, offset_z, offset_c);
        offset_c += MIDDLE_KNOT_CONSTRAINTS;
    }

    offset_z = KNOT_SIZE*(N_KNOTS-1);
    Load_End_Knot_Columns(offset_z, offset_c);
}

void Parse_Checklist_S(){
    int k_ic, k_nz, k_tot, k_j;
    int *ptr;
    int i, j;
    int P_temp[SYSTEM_SIZE];

    // Reset the Checklist so it looks like it started before being permuted
    memcpy(P_temp, P, SYSTEM_SIZE*sizeof(int));

    for (i = 0; i < SYSTEM_SIZE; i++){
        while(P_temp[i] != i){
            j = P_temp[i];
            P_temp[i] = P_temp[j];
            P_temp[j] = j;

            ptr = Checklist_S[i];
            Checklist_S[i] = Checklist_S[j];
            Checklist_S[j] = ptr;
        }
    }

    A_S_ids = malloc(SYSTEM_SIZE*sizeof(int*));
    Zero_cols = malloc(SYSTEM_SIZE*sizeof(int*));
    J_S_ids = malloc(N_DECISION_VARIABLES*sizeof(int*));
    Lambda_ids = malloc(N_DECISION_VARIABLES*sizeof(int*));

    for (i = 0; i < SYSTEM_SIZE; i++){

        for (j = 0; j < SYSTEM_SIZE; j++){
            if (Checklist_S[i][j] != 0)
                Sparse_A_Size[i]++;

            if (Checklist_S[i][j] == 1){
                Col_ic_S[i]++;
                if (i < N_DECISION_VARIABLES && j >= N_DECISION_VARIABLES)
                    Col_nz_J[i]++;
            }

            else if (Checklist_S[i][j] == 2)
                Col_nz_S[i]++;
        }
        A_S_ids[i] = malloc(Sparse_A_Size[i]*sizeof(int));
        Zero_cols[i] = malloc(Col_nz_S[i]*sizeof(int));
        if (i < N_DECISION_VARIABLES){
            J_S_ids[i] = malloc(Col_nz_J[i]*sizeof(int));
            Lambda_ids[i] = malloc(Col_nz_J[i]*sizeof(int));
        }

        k_ic = 0;
        k_nz = 0;
        k_tot = 0;
        k_j = 0;
        for (j = 0; j < SYSTEM_SIZE; j++){

            if (Checklist_S[i][j] == 1){
                k_ic++;
                if (i < N_DECISION_VARIABLES && j >= N_DECISION_VARIABLES){
                    Lambda_ids[i][k_j] = j - N_DECISION_VARIABLES;
                    J_S_ids[i][k_j++] = k_tot;
                }
            }

            else if (Checklist_S[i][j] == 2)
                Zero_cols[i][k_nz++] = k_tot;

            if (Checklist_S[i][j] != 0)
                A_S_ids[i][k_tot++] = j;
        }

    }

    Get_IC_Cols();
}

void Fill_A_S_Zeros(){
    int i, j;
    for (i = 0; i < SYSTEM_SIZE; i++)
        for (j = 0; j < Col_nz_S[i]; j++)
            A_S[i][Zero_cols[i][j]] = 0.0;
}

void Convert_id_U_S(){
    int i, j, k;
    for (i = 0; i < SYSTEM_SIZE; i++){
        for (j = 0; j < Col_nz_U[i]; j++)
            Col_ids_S1[i][j] = Search_For_Sparse_Column(P[i], Col_ids_U[i][j]);
        for (j = 0; j < Col_nz_L[i]; j++)
            Col_ids_S2[i][j] = Search_For_Sparse_Column(P[i], Col_ids_L[i][j]);
        for (j = 0; j < Row_nz[i]; j++)
            for (k = 0; k < Col_nz_U[i]; k++)
                Col_ids_S3[i][j][k] = Search_For_Sparse_Column(Global_Row[i][j], Col_ids_S3[i][j][k]);
    }

}

double *Get_z(){
    return &z[0];
}

double *Get_t(){
    return &t[0];
}