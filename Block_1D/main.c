#include <stdint.h>
#include <stdlib.h>
#include <stdio.h>
#include <math.h>
#include <time.h>
#include <string.h>
#include "LU.h"

void Run_SQP();
void Load_Problem_Data();
void Load_Jacobian();
void Load_Hessian();
void Load_Equalities();
void Load_Gradient();
double Select_Alpha(double *delta_z);
double Get_Cost(double *z);
void Reset_L();

#define N_STATES 2
#define N_INPUTS 1
#define N_KNOTS 300
#define N_LBS 1
#define N_UBS 1
#define FIRST_KNOT_CONSTRAINTS 6
#define MIDDLE_KNOT_CONSTRAINTS 4
#define END_KNOT_CONSTRAINTS 2

#define KNOT_SIZE (N_STATES + N_INPUTS + N_LBS + N_UBS)
#define N_DECISION_VARIABLES (N_KNOTS * KNOT_SIZE)
#define N_CONSTRAINTS (FIRST_KNOT_CONSTRAINTS + MIDDLE_KNOT_CONSTRAINTS*(N_KNOTS-2) + END_KNOT_CONSTRAINTS)
#define SYSTEM_SIZE (N_DECISION_VARIABLES + N_CONSTRAINTS)

int plot_flag = 1;
double ic[2] = {0.0, 0.0};
double xd[2] = {2.0, 0.4};
double lb = -2.0;
double ub = 2.0;
double tf = 2.0;
double mu = 1.0e-7;

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

int Row_nz[SYSTEM_SIZE];
int Col_nz_U[SYSTEM_SIZE];
int Col_nz_L[SYSTEM_SIZE];

int **Row_ids;
// These are the column indices for each row above the diagonal (elimination, backwards substitution)
int **Col_ids_U;
// These are the column indices for each row below the diagonal (forwards substitution)
int **Col_ids_L;
// Lower factorization of A
double **A_L;
// Upper factorization of A, includes diagonal
double **A_U;

typedef struct{
    double dx[N_KNOTS][N_STATES];
}params;
params knot_params = {0};
struct timespec time_load, time_elim, time_subs, time_extract, time_alpha, time_update;

int main(){
    int i, j;
    struct timespec time_start, time_end;

    // Fill t and h
    double dt = tf / (double)(N_KNOTS - 1);
    for (i = 1; i < N_KNOTS; i++){
        t[i] = t[i-1] + dt;
        h[i-1] = dt; 
    }

    // Set initial conditions
    z[0] = ic[0];
    z[1] = ic[1];

    // Set slack variables
    for (i = 0; i < N_KNOTS; i++){
        j = KNOT_SIZE*i;
        z[j+3] = z[j+2] - lb;
        z[j+4] = ub - z[j+2];
    }

    // Allocate A
    A = malloc(SYSTEM_SIZE*sizeof(double*));
    if (A == NULL) return 0;
    for (i = 0; i < SYSTEM_SIZE; i++){
        A[i] = calloc(SYSTEM_SIZE, sizeof(double));
        if (A[i] == NULL) return 0;
    }

    // Find Pivots and Extents
    Load_Problem_Data();
    Get_Pivots_Extents(A, SYSTEM_SIZE, 1.0e-20, P, P_max, Row_nz, Col_nz_U, Col_nz_L);
    for (i = 0; i < SYSTEM_SIZE; i++){
        memset(A[i], 0, SYSTEM_SIZE*sizeof(double));
    }

    // Allocate stop indices
    Row_ids = malloc(SYSTEM_SIZE*sizeof(int*));
    if (Row_ids == NULL) return 0;
    for (i = 0; i < SYSTEM_SIZE; i++){
        Row_ids[i] = malloc(Row_nz[i]*sizeof(int));
        if (Row_ids[i] == NULL) return 0;
    }

    Col_ids_U = malloc(SYSTEM_SIZE*sizeof(int*));
    A_U = malloc(SYSTEM_SIZE*sizeof(double*));
    for (i = 0; i < SYSTEM_SIZE; i++){
        Col_ids_U[i] = malloc(Col_nz_U[i]*sizeof(int));
        A_U[i] = malloc((Col_nz_U[i]+1)*sizeof(double));
    }

    Col_ids_L = malloc(SYSTEM_SIZE*sizeof(int*));
    A_L = malloc(SYSTEM_SIZE*sizeof(double*));
    for (i = 0; i < SYSTEM_SIZE; i++){
        Col_ids_L[i] = malloc(Col_nz_L[i]*sizeof(int));
        A_L[P[i]] = malloc(Col_nz_L[i]*sizeof(double));
    }

    // Find stop indices
    Load_Problem_Data();
    Get_Stops(A, SYSTEM_SIZE, P_max, Row_ids, Col_ids_U, Col_ids_L);
    for (i = 0; i < SYSTEM_SIZE; i++){
        memset(A[i], 0, SYSTEM_SIZE*sizeof(double));
    }

    timespec_get(&time_start, TIME_UTC);

    // Run problem
    Run_SQP();

    timespec_get(&time_end, TIME_UTC);

    // Free A
    for (i = 0; i < SYSTEM_SIZE; i++){
        free(A[i]);
        free(Row_ids[i]);
        free(Col_ids_U[i]);
        free(Col_ids_L[i]);
        free(A_L[i]);
        free(A_U[i]);
    }
    free(A);
    free(Row_ids);
    free(Col_ids_U);
    free(Col_ids_L);
    free(A_L);
    free(A_U);

    // Log data
    if (plot_flag){
        FILE *log = fopen("Block1D_log.txt", "w");
        for (i = 0; i < N_KNOTS; i++){
            int offset_z = i*KNOT_SIZE;
            fprintf(log, "%6.3f  %6.3f  %6.3f  %6.3f\n", t[i], z[offset_z], z[offset_z+1], z[offset_z+2]);
        }
        fclose(log);
    }
    // Plot results
    if (plot_flag) system("gnuplot plotter.plt");

    double time_elapsed = (double)(time_end.tv_sec - time_start.tv_sec) + (double)(time_end.tv_nsec - time_start.tv_nsec)/1.0e9;
    double load_time = (double)(time_load.tv_sec - time_start.tv_sec) + (double)(time_load.tv_nsec - time_start.tv_nsec)/1.0e9;
    double elim_time = (double)(time_elim.tv_sec - time_load.tv_sec) + (double)(time_elim.tv_nsec - time_load.tv_nsec)/1.0e9;
    double subs_time = (double)(time_subs.tv_sec - time_elim.tv_sec) + (double)(time_subs.tv_nsec - time_elim.tv_nsec)/1.0e9;
    double extr_time = (double)(time_extract.tv_sec - time_subs.tv_sec) + (double)(time_extract.tv_nsec - time_subs.tv_nsec)/1.0e9;
    double alpha_time = (double)(time_alpha.tv_sec - time_extract.tv_sec) + (double)(time_alpha.tv_nsec - time_extract.tv_nsec)/1.0e9;
    double upd_time = (double)(time_update.tv_sec - time_alpha.tv_sec) + (double)(time_update.tv_nsec - time_alpha.tv_nsec)/1.0e9;
    printf("%.9f seconds elapsed\n", time_elapsed);
    printf("%.9f seconds loading\n", load_time);
    printf("%.9f seconds eliminating\n", elim_time);
    printf("%.9f seconds substituting\n", subs_time);
    printf("%.9f seconds extracting\n", extr_time);
    printf("%.9f seconds alpha\n", alpha_time);
    printf("%.9f seconds updating\n", upd_time);
    return 1;
}

void Run_SQP(){
    double delta_z[N_DECISION_VARIABLES];
    double lambda_new[N_CONSTRAINTS];
    double delta_lambda[N_CONSTRAINTS];
    const int max_iterations = 2;
    int i, iterations;
    
    for (iterations = 0; iterations < max_iterations; iterations++){
        // Load A and b matrices
        // for (i = 0; i < SYSTEM_SIZE; i++)
            // memset(A[i], 0, SYSTEM_SIZE*sizeof(double));
        Load_Problem_Data();
        timespec_get(&time_load, TIME_UTC);

        // Setup A in form A=(L-E)+U by partial pivoting and eliminating
        (void)LUDecompose(A, SYSTEM_SIZE, P_max,
                            Row_nz, Row_ids,
                            Col_nz_U, Col_ids_U,
                            A_L, A_U);
        timespec_get(&time_elim, TIME_UTC);

        // Perform forward and backward substitution to solve for x
        LUPSolve(A, P, b, SYSTEM_SIZE, x,
                Col_nz_U, Col_ids_U,
                Col_nz_L, Col_ids_L,
                A_L, A_U);
        // Need to swap back A_L after solving so it is ready for next iteration
        Reset_L();
        timespec_get(&time_subs, TIME_UTC);

        // Extract change in decision variables and new lagrange multipliers from x
        memcpy(delta_z, x, N_DECISION_VARIABLES*sizeof(double));
        memcpy(lambda_new, &x[N_DECISION_VARIABLES], N_CONSTRAINTS*sizeof(double));
        for (i = 0; i < N_CONSTRAINTS; i++)
            delta_lambda[i] = lambda_new[i] - lambda[i];
        timespec_get(&time_extract, TIME_UTC);

        // Check slack variables and perform ternary search to determine what value of alpha to use
        double alpha = Select_Alpha(delta_z);
        timespec_get(&time_alpha, TIME_UTC);

        // Update decision variables and lagrange multipliers with selected alpha
        for (i = 0; i < N_DECISION_VARIABLES; i++)
            z[i] += alpha*delta_z[i];
        
        for (i = 0; i < N_CONSTRAINTS; i++)
            lambda[i] += alpha*delta_lambda[i];
        timespec_get(&time_update, TIME_UTC);

    }
}

void Reset_L(){
    double *ptr;
    int i, j;
    int P_temp[SYSTEM_SIZE];
    memcpy(P_temp, P, SYSTEM_SIZE*sizeof(int));

    for (i = 0; i < SYSTEM_SIZE; i++){
        while(P_temp[i] != i){
            j = P_temp[i];
            P_temp[i] = P_temp[j];
            P_temp[j] = j;

            ptr = A_L[i];
            A_L[i] = A_L[j];
            A_L[j] = ptr;
        }
    }
}

double First_Knot_Cost(double *z_knot){
    double contribution = -mu*(log(z_knot[3]) + log(z_knot[4]));
    return contribution;
}

double Middle_Knot_Cost(double *z_knot){
    double contribution = -mu*(log(z_knot[3]) + log(z_knot[4]));
    return contribution;
}

double End_Knot_Cost(double *z_knot){
    double contribution = pow(z_knot[0] - xd[0], 2) + pow(z_knot[1] - xd[1], 2) - mu*(log(z_knot[3]) + log(z_knot[4]));
    return contribution;
}

double Get_Cost(double *z){
    int i, offset_z;
    double cost = 0;

    offset_z = 0;
    cost += First_Knot_Cost(&z[offset_z]);

    for (i = 1; i < N_KNOTS - 1; i++){
        offset_z = i*KNOT_SIZE;
        cost += Middle_Knot_Cost(&z[offset_z]);
    }

    offset_z = KNOT_SIZE*(N_KNOTS-1);
    cost += End_Knot_Cost(&z[offset_z]);

    return cost;
}

double Select_Alpha(double *delta_z){
    double z_new[N_DECISION_VARIABLES];
    double alpha_check;
    double alpha = 1.0;
    const double tau = 0.995;
    int i, offset_z;
    // First check the maximum alpha that will keep our slack variables positive
    for (i = 0; i < N_DECISION_VARIABLES; i++)
        z_new[i] = z[i] + delta_z[i];

    for (i = 0; i < N_KNOTS; i++){
        offset_z = i*KNOT_SIZE;
        if (z_new[offset_z+3] <= 0.0 && delta_z[offset_z+3] < 0.0){
            alpha_check = - z[offset_z+3] / delta_z[offset_z+3];
            if (alpha_check < alpha)
                alpha = alpha_check;
        }
        if (z_new[offset_z+4] <= 0.0 && delta_z[offset_z+4] < 0.0){
            alpha_check = - z[offset_z+4] / delta_z[offset_z+4];
            if (alpha_check < alpha)
                alpha = alpha_check;
        }

    }

    alpha *= tau;

    double alow = 0.0;
    double ahigh = alpha;
    double amid1, amid2, cost1, cost2;
    double zmid1[N_DECISION_VARIABLES];
    double zmid2[N_DECISION_VARIABLES];
    int j;
    const int max_ternary = 12;
    for (i = 0; i < max_ternary; i++){
        amid1 = alow + (1.0/3.0)*(ahigh - alow);
        amid2 = alow + (2.0/3.0)*(ahigh - alow);

        for (j = 0; j < N_DECISION_VARIABLES; j++){
            zmid1[j] = z[j] + amid1*delta_z[j];
            zmid2[j] = z[j] + amid2*delta_z[j];
        }
        
        cost1 = Get_Cost(zmid1);
        cost2 = Get_Cost(zmid2);

        if (cost1 <= cost2)
            ahigh = amid2;
        else
            alow = amid1;
    }

    alpha = (amid1 + amid2)/2.0;

    return alpha;
}

void Load_Problem_Data(){
    // Update and load the gradient of the cost function wrt the decision variables
    Load_Gradient();

    // Evaluate and load the equality constraints
    // This also updates the knot parameters
    Load_Equalities();

    // Evaluate and load the jacobian of the constraints wrt the decision variables
    Load_Jacobian();
    
    // Evaluate and load the hessian of the lagrangian wrt the decision variables
    Load_Hessian();
}

void First_Knot_Jacobian(int knot, int offset_z, int offset_c){
    int i = N_DECISION_VARIABLES;
    int j = i + offset_c;
    int k = offset_z;
    const double c1 = -h[knot] / 2.0;

    A[j][k] = -1.0;
    A[j][k+1] = c1;
    A[j][k+5] = 1.0;
    A[j][k+6] = c1;
    A[k][j] = A[j][k];
    A[k+1][j] = A[j][k+1];
    A[k+5][j] = A[j][k+5];
    A[k+6][j] = A[j][k+6];

    A[j+1][k+1] = -1.0;
    A[j+1][k+2] = c1;
    A[j+1][k+6] = 1.0;
    A[j+1][k+7] = c1;
    A[k+1][j+1] = A[j+1][k+1];
    A[k+2][j+1] = A[j+1][k+2];
    A[k+6][j+1] = A[j+1][k+6];
    A[k+7][j+1] = A[j+1][k+7];

    A[j+2][k+2] = -1.0;
    A[j+2][k+3] = 1.0;
    A[k+2][j+2] = A[j+2][k+2];
    A[k+3][j+2] = A[j+2][k+3];

    A[j+3][k+2] = 1.0;
    A[j+3][k+4] = 1.0;
    A[k+2][j+3] = A[j+3][k+2];
    A[k+4][j+3] = A[j+3][k+4];

    A[j+4][k] = -1.0;
    A[k][j+4] = A[j+4][k];

    A[j+5][k+1] = -1.0;
    A[k+1][j+5] = A[j+5][k+1];
}

void Middle_Knot_Jacobian(int knot, int offset_z, int offset_c){
    int i = N_DECISION_VARIABLES;
    int j = i + offset_c;
    int k = offset_z;
    const double c1 = -h[knot] / 2.0;

    A[j][k] = -1.0;
    A[j][k+1] = c1;
    A[j][k+5] = 1.0;
    A[j][k+6] = c1;
    A[k][j] = A[j][k];
    A[k+1][j] = A[j][k+1];
    A[k+5][j] = A[j][k+5];
    A[k+6][j] = A[j][k+6];

    A[j+1][k+1] = -1.0;
    A[j+1][k+2] = c1;
    A[j+1][k+6] = 1.0;
    A[j+1][k+7] = c1;
    A[k+1][j+1] = A[j+1][k+1];
    A[k+2][j+1] = A[j+1][k+2];
    A[k+6][j+1] = A[j+1][k+6];
    A[k+7][j+1] = A[j+1][k+7];

    A[j+2][k+2] = -1.0;
    A[j+2][k+3] = 1.0;
    A[k+2][j+2] = A[j+2][k+2];
    A[k+3][j+2] = A[j+2][k+3];

    A[j+3][k+2] = 1.0;
    A[j+3][k+4] = 1.0;
    A[k+2][j+3] = A[j+3][k+2];
    A[k+4][j+3] = A[j+3][k+4];

}

void End_Knot_Jacobian(int knot, int offset_z, int offset_c){
    int i = N_DECISION_VARIABLES;
    int j = i + offset_c;
    int k = offset_z;

    A[j][k+2] = -1.0;
    A[j][k+3] = 1.0;
    A[k+2][j] = A[j][k+2];
    A[k+3][j] = A[j][k+3];

    A[j+1][k+2] = 1.0;
    A[j+1][k+4] = 1.0;
    A[k+2][j+1] = A[j+1][k+2];
    A[k+4][j+1] = A[j+1][k+4];

}

void Load_Jacobian(){
    int i, offset_z, offset_c;
    offset_z = 0;
    offset_c = 0;
    First_Knot_Jacobian(0, offset_z, offset_c);
    offset_c += FIRST_KNOT_CONSTRAINTS;

    for (i = 1; i < N_KNOTS - 1; i++){
        offset_z = i*KNOT_SIZE;
        Middle_Knot_Jacobian(i, offset_z, offset_c);
        offset_c += MIDDLE_KNOT_CONSTRAINTS;
    }

    offset_z = KNOT_SIZE*(N_KNOTS-1);
    End_Knot_Jacobian(N_KNOTS-1, offset_z, offset_c);
}

void First_Knot_Hessian(int knot, int offset_z){
    int i = offset_z;
    
    A[i+3][i+3] = mu / pow(z[i+3], 2);
    A[i+4][i+4] = mu / pow(z[i+4], 2);
}

void Middle_Knot_Hessian(int knot, int offset_z){
    int i = offset_z;
    
    A[i+3][i+3] = mu / pow(z[i+3], 2);
    A[i+4][i+4] = mu / pow(z[i+4], 2);
}

void End_Knot_Hessian(int knot, int offset_z){
    int i = offset_z;
    
    A[i][i] = 2.0;
    A[i+1][i+1] = 2.0;
    A[i+3][i+3] = mu / pow(z[i+3], 2);
    A[i+4][i+4] = mu / pow(z[i+4], 2);
}

void Load_Hessian(){
    int i, offset_z;
    offset_z = 0;
    First_Knot_Hessian(0, offset_z);

    for (i = 1; i < N_KNOTS - 1; i++){
        offset_z = i*KNOT_SIZE;
        Middle_Knot_Hessian(i, offset_z);
    }

    offset_z = KNOT_SIZE*(N_KNOTS-1);
    End_Knot_Hessian(N_KNOTS-1, offset_z);
}

void Get_Knot_Params(int offset, int knot){
    double v = z[offset+1];
    double a = z[offset+2];
    
    knot_params.dx[knot][0] = v;
    knot_params.dx[knot][1] = a;
}

void First_Knot_Eq(int knot, int offset_z, int offset_c){
    // p1 = z[offset_z];
    // v1 = z[offset_z+1];
    // a1 = z[offset_z+2];
    // s_a_lb = z[offset_z+3];
    // s_a_ub = z[offset_z+4];
    // p2 = z[offset_z+5];
    // v2 = z[offset_z+6];

    // Defects
    b[N_DECISION_VARIABLES+offset_c] = 
        0.5*h[knot]*(knot_params.dx[knot][0] + knot_params.dx[knot+1][0]) + z[offset_z] - z[offset_z+5];
    b[N_DECISION_VARIABLES+offset_c+1] = 
        0.5*h[knot]*(knot_params.dx[knot][1] + knot_params.dx[knot+1][1]) + z[offset_z+1] - z[offset_z+6];
    // Lower bounds
    b[N_DECISION_VARIABLES+offset_c+2] = z[offset_z+2] - lb - z[offset_z+3];
    // Upper bounds
    b[N_DECISION_VARIABLES+offset_c+3] = ub - z[offset_z+2] - z[offset_z+4];
    // Initial conditions
    b[N_DECISION_VARIABLES+offset_c+4] = z[offset_z] - ic[0];
    b[N_DECISION_VARIABLES+offset_c+5] = z[offset_z+1] - ic[1];
}

void Middle_Knot_Eq(int knot, int offset_z, int offset_c){
    // p1 = z[offset_z];
    // v1 = z[offset_z+1];
    // a1 = z[offset_z+2];
    // s_a_lb = z[offset_z+3];
    // s_a_ub = z[offset_z+4];
    // p2 = z[offset_z+5];
    // v2 = z[offset_z+6];

    // Defects
    b[N_DECISION_VARIABLES+offset_c] = 
        0.5*h[knot]*(knot_params.dx[knot][0] + knot_params.dx[knot+1][0]) + z[offset_z] - z[offset_z+5];
    b[N_DECISION_VARIABLES+offset_c+1] = 
        0.5*h[knot]*(knot_params.dx[knot][1] + knot_params.dx[knot+1][1]) + z[offset_z+1] - z[offset_z+6];
    // Lower bounds
    b[N_DECISION_VARIABLES+offset_c+2] = z[offset_z+2] - lb - z[offset_z+3];
    // Upper bounds
    b[N_DECISION_VARIABLES+offset_c+3] = ub - z[offset_z+2] - z[offset_z+4];

}

void End_Knot_Eq(int knot, int offset_z, int offset_c){
    // p1 = z[offset_z];
    // v1 = z[offset_z+1];
    // a1 = z[offset_z+2];
    // s_a_lb = z[offset_z+3];
    // s_a_ub = z[offset_z+4];

    // Lower bounds
    b[N_DECISION_VARIABLES+offset_c] = z[offset_z+2] - lb - z[offset_z+3];
    // Upper bounds
    b[N_DECISION_VARIABLES+offset_c+1] = ub - z[offset_z+2] - z[offset_z+4];
}

void Load_Equalities(){
    int i, offset_z, offset_c;

    for (i = 0; i < N_KNOTS; i++){
        offset_z = KNOT_SIZE*i;
        Get_Knot_Params(offset_z, i);
    }

    offset_z = 0;
    offset_c = 0;
    First_Knot_Eq(0, offset_z, offset_c);
    offset_c += FIRST_KNOT_CONSTRAINTS;

    for (i = 1; i < N_KNOTS - 1; i++){
        offset_z = i*KNOT_SIZE;
        Middle_Knot_Eq(i, offset_z, offset_c);
        offset_c += MIDDLE_KNOT_CONSTRAINTS;
    }

    offset_z = KNOT_SIZE*(N_KNOTS-1);
    End_Knot_Eq(N_KNOTS-1, offset_z, offset_c);
}

void First_Knot_Grad(int offset){
    b[offset+3] = mu / z[offset+3];
    b[offset+4] = mu / z[offset+4];
}

void Middle_Knot_Grad(int offset){
    b[offset+3] = mu / z[offset+3];
    b[offset+4] = mu / z[offset+4];
}

void End_Knot_Grad(int offset){
    b[offset] = 2.0*(xd[0] - z[offset]);
    b[offset+1] = 2.0*(xd[1] - z[offset+1]);
    b[offset+3] = mu / z[offset+3];
    b[offset+4] = mu / z[offset+4];
}

void Load_Gradient(){
    int i, offset;

    offset = 0;
    First_Knot_Grad(offset);

    for (i = 1; i < N_KNOTS - 1; i++){
        offset = KNOT_SIZE*i;
        Middle_Knot_Grad(offset);
    }

    offset = KNOT_SIZE*(N_KNOTS-1);
    End_Knot_Grad(offset);
}
