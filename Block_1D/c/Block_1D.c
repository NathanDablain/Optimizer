#include <math.h>
#include <stdlib.h>
#include <stdio.h>
#include "Block_1D.h"
#include "SQP.h"

double ic[2] = {0.0, 0.0};
double xd[2] = {2.0, 0.4};
double lb = -2.0;
double ub = 2.0;
double tf = 2.0;
double mu = 1.0e-7;

typedef struct{
    double dx[N_KNOTS][N_STATES];
}params;

typedef struct{
    int k1[FIRST_KNOT_CONTRIBUTION];
    int km[N_KNOTS-2][MIDDLE_KNOT_CONTRIBUTION];
    int ke[END_KNOT_CONTRIBUTION];
    params knot_params;
}knots;

knots problem;

double Get_tf(){
    return tf;
}

void Load_ic(double *z){
    z[0] = ic[0];
    z[1] = ic[1];
}

void Load_Problem_Data(double **A, double *b, double *z, double *h, bool sparse){
    // Update and load the gradient of the cost function wrt the decision variables
    Load_Gradient(b, z);

    // Evaluate and load the equality constraints
    // This also updates the knot parameters
    Load_Equalities(b, z, h);

    // Evaluate and load the jacobian of the constraints wrt the decision variables
    Load_Jacobian(A, h, sparse);
    
    // Evaluate and load the hessian of the lagrangian wrt the decision variables
    Load_Hessian(A, z, sparse);

}

double First_Knot_Cost(double *z_knot){
    double contribution = -mu*(log(z_knot[2] - lb) +
                               log(ub - z_knot[2]));
    return contribution;
}

double Middle_Knot_Cost(double *z_knot){
    double contribution = -mu*(log(z_knot[2] - lb) +
                               log(ub - z_knot[2]));
    return contribution;
}

double End_Knot_Cost(double *z_knot){
    double contribution = pow(z_knot[0] - xd[0], 2) +
                          pow(z_knot[1] - xd[1], 2) -
                          mu*(log(z_knot[2] - lb) +
                              log(ub - z_knot[2]));
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

double Select_Alpha(double *z, double *delta_z){
    double z_new[N_DECISION_VARIABLES];
    double slack_new[N_LBS + N_UBS];
    double slack_cur[N_LBS + N_UBS];
    double alpha_check;
    double alpha = 1.0;
    const double tau = 0.995;
    int i, j, offset_z;
    // First check the maximum alpha that will maintain positive slack
    for (i = 0; i < N_DECISION_VARIABLES; i++)
        z_new[i] = z[i] + delta_z[i];

    for (i = 0; i < N_KNOTS; i++){
        offset_z = i*KNOT_SIZE;
        slack_new[0] = z_new[offset_z+2] - lb;
        slack_new[1] = ub - z_new[offset_z+2];
        slack_cur[0] = z[offset_z+2] - lb;
        slack_cur[1] = ub - z[offset_z+2];
        for (j = 0; j < N_LBS + N_UBS; j++){
            if (slack_new[j] <= 0.0) {
                alpha_check = -slack_cur[j] / (slack_new[j] - slack_cur[j]);
                if (alpha_check < alpha)
                    alpha = alpha_check;
            }
        }

    }

    alpha *= tau;

    double alow = 0.0;
    double ahigh = alpha;
    double amid1, amid2, cost1, cost2;
    double zmid1[N_DECISION_VARIABLES];
    double zmid2[N_DECISION_VARIABLES];
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

void First_Knot_Jacobian(double **A, double *h, int knot, int offset_z, int offset_c, bool sparse){
    int i = N_DECISION_VARIABLES;
    int j = i + offset_c;
    int k = offset_z;
    const double c1 = -h[knot] / 2.0;

    if (sparse){
        A[j][problem.k1[0]] = -1.0;
        A[j][problem.k1[1]] = c1;
        A[j][problem.k1[2]] = 1.0;
        A[j][problem.k1[3]] = c1;
        A[k][problem.k1[4]] = A[j][problem.k1[0]];
        A[k+1][problem.k1[5]] = A[j][problem.k1[1]];
        A[k+3][problem.k1[6]] = A[j][problem.k1[2]];
        A[k+4][problem.k1[7]] = A[j][problem.k1[3]];

        A[j+1][problem.k1[8]] = -1.0;
        A[j+1][problem.k1[9]] = c1;
        A[j+1][problem.k1[10]] = 1.0;
        A[j+1][problem.k1[11]] = c1;
        A[k+1][problem.k1[12]] = A[j+1][problem.k1[8]];
        A[k+2][problem.k1[13]] = A[j+1][problem.k1[9]];
        A[k+4][problem.k1[14]] = A[j+1][problem.k1[10]];
        A[k+5][problem.k1[15]] = A[j+1][problem.k1[11]];

        A[j+2][problem.k1[16]] = -1.0;
        A[k][problem.k1[17]] = A[j+2][problem.k1[16]];

        A[j+3][problem.k1[18]] = -1.0;
        A[k+1][problem.k1[19]] = A[j+3][problem.k1[18]];
    }
    else{
        A[j][k] = -1.0;
        A[j][k+1] = c1;
        A[j][k+3] = 1.0;
        A[j][k+4] = c1;
        A[k][j] = A[j][k];
        A[k+1][j] = A[j][k+1];
        A[k+3][j] = A[j][k+3];
        A[k+4][j] = A[j][k+4];

        A[j+1][k+1] = -1.0;
        A[j+1][k+2] = c1;
        A[j+1][k+4] = 1.0;
        A[j+1][k+5] = c1;
        A[k+1][j+1] = A[j+1][k+1];
        A[k+2][j+1] = A[j+1][k+2];
        A[k+4][j+1] = A[j+1][k+4];
        A[k+5][j+1] = A[j+1][k+5];

        A[j+2][k] = -1.0;
        A[k][j+2] = A[j+2][k];

        A[j+3][k+1] = -1.0;
        A[k+1][j+3] = A[j+3][k+1];
    }
}

void Middle_Knot_Jacobian(double **A, double *h, int knot, int offset_z, int offset_c, bool sparse){
    int i = N_DECISION_VARIABLES;
    int j = i + offset_c;
    int k = offset_z;
    const double c1 = -h[knot] / 2.0;

    if (sparse){
        A[j][problem.km[knot-1][0]] = -1.0;
        A[j][problem.km[knot-1][1]] = c1;
        A[j][problem.km[knot-1][2]] = 1.0;
        A[j][problem.km[knot-1][3]] = c1;
        A[k][problem.km[knot-1][4]] = A[j][problem.km[knot-1][0]];
        A[k+1][problem.km[knot-1][5]] = A[j][problem.km[knot-1][1]];
        A[k+3][problem.km[knot-1][6]] = A[j][problem.km[knot-1][2]];
        A[k+4][problem.km[knot-1][7]] = A[j][problem.km[knot-1][3]];

        A[j+1][problem.km[knot-1][8]] = -1.0;
        A[j+1][problem.km[knot-1][9]] = c1;
        A[j+1][problem.km[knot-1][10]] = 1.0;
        A[j+1][problem.km[knot-1][11]] = c1;
        A[k+1][problem.km[knot-1][12]] = A[j+1][problem.km[knot-1][8]];
        A[k+2][problem.km[knot-1][13]] = A[j+1][problem.km[knot-1][9]];
        A[k+4][problem.km[knot-1][14]] = A[j+1][problem.km[knot-1][10]];
        A[k+5][problem.km[knot-1][15]] = A[j+1][problem.km[knot-1][11]];

    }
    else{
        A[j][k] = -1.0;
        A[j][k+1] = c1;
        A[j][k+3] = 1.0;
        A[j][k+4] = c1;
        A[k][j] = A[j][k];
        A[k+1][j] = A[j][k+1];
        A[k+3][j] = A[j][k+3];
        A[k+4][j] = A[j][k+4];

        A[j+1][k+1] = -1.0;
        A[j+1][k+2] = c1;
        A[j+1][k+4] = 1.0;
        A[j+1][k+5] = c1;
        A[k+1][j+1] = A[j+1][k+1];
        A[k+2][j+1] = A[j+1][k+2];
        A[k+4][j+1] = A[j+1][k+4];
        A[k+5][j+1] = A[j+1][k+5];

    }
}

void End_Knot_Jacobian(double **A, int knot, int offset_z, int offset_c, bool sparse){

}

void Load_Jacobian(double **A, double *h, bool sparse){
    int i, offset_z, offset_c;
    offset_z = 0;
    offset_c = 0;

    First_Knot_Jacobian(A, h, 0, offset_z, offset_c, sparse);
    offset_c += FIRST_KNOT_CONSTRAINTS;

    for (i = 1; i < N_KNOTS - 1; i++){
        offset_z = i*KNOT_SIZE;
        Middle_Knot_Jacobian(A, h, i, offset_z, offset_c, sparse);
        offset_c += MIDDLE_KNOT_CONSTRAINTS;
    }

    offset_z = KNOT_SIZE*(N_KNOTS-1);
    End_Knot_Jacobian(A, N_KNOTS-1, offset_z, offset_c, sparse);
}

void First_Knot_Hessian(double **A, double *z, int knot, int offset_z, bool sparse){
    int i = offset_z;
    
    if (sparse){
        A[i+2][problem.k1[20]] = (mu/pow(z[i+2] - ub, 2)) + (mu/pow(lb - z[i+2], 2));
    }
    else{
        A[i+2][i+2] = (mu/pow(z[i+2] - ub, 2)) + (mu/pow(lb - z[i+2], 2));
    }
}

void Middle_Knot_Hessian(double **A, double *z, int knot, int offset_z, bool sparse){
    int i = offset_z;
    if (sparse){
        A[i+2][problem.km[knot-1][16]] = (mu/pow(z[i+2] - ub, 2)) + (mu/pow(lb - z[i+2], 2));
    }
    else{
        A[i+2][i+2] = (mu/pow(z[i+2] - ub, 2)) + (mu/pow(lb - z[i+2], 2));
    }
}

void End_Knot_Hessian(double **A, double *z, int knot, int offset_z, bool sparse){
    int i = offset_z;
    
    if (sparse){
        A[i][problem.ke[0]] = 2.0;
        A[i+1][problem.ke[1]] = 2.0;
        A[i+2][problem.ke[2]] = (mu/pow(z[i+2] - ub, 2)) + (mu/pow(lb - z[i+2], 2));
    }
    else{
        A[i][i] = 2.0;
        A[i+1][i+1] = 2.0;
        A[i+2][i+2] = (mu/pow(z[i+2] - ub, 2)) + (mu/pow(lb - z[i+2], 2));
    }

}

void Load_Hessian(double **A, double *z, bool sparse){
    int i, offset_z;
    offset_z = 0;

    First_Knot_Hessian(A, z, 0, offset_z, sparse);

    for (i = 1; i < N_KNOTS - 1; i++){
        offset_z = i*KNOT_SIZE;
        Middle_Knot_Hessian(A, z, i, offset_z, sparse);
    }

    offset_z = KNOT_SIZE*(N_KNOTS-1);
    End_Knot_Hessian(A, z, N_KNOTS-1, offset_z, sparse);
}

void Get_Knot_Params(double *z, int offset, int knot){
    double v = z[offset+1];
    double a = z[offset+2];
    
    problem.knot_params.dx[knot][0] = v;
    problem.knot_params.dx[knot][1] = a;
}

void First_Knot_Eq(double *b, double *z, double *h, int knot, int offset_z, int offset_c){
    // p1 = z[offset_z];
    // v1 = z[offset_z+1];
    // a1 = z[offset_z+2];
    // p2 = z[offset_z+5];
    // v2 = z[offset_z+6];

    // Defects
    b[N_DECISION_VARIABLES+offset_c] = 
        0.5*h[knot]*(problem.knot_params.dx[knot][0] + problem.knot_params.dx[knot+1][0]) + z[offset_z] - z[offset_z+3];
    b[N_DECISION_VARIABLES+offset_c+1] = 
        0.5*h[knot]*(problem.knot_params.dx[knot][1] + problem.knot_params.dx[knot+1][1]) + z[offset_z+1] - z[offset_z+4];
    // Initial conditions
    b[N_DECISION_VARIABLES+offset_c+2] = z[offset_z] - ic[0];
    b[N_DECISION_VARIABLES+offset_c+3] = z[offset_z+1] - ic[1];
}

void Middle_Knot_Eq(double *b, double *z, double *h, int knot, int offset_z, int offset_c){
    // p1 = z[offset_z];
    // v1 = z[offset_z+1];
    // a1 = z[offset_z+2];
    // p2 = z[offset_z+3];
    // v2 = z[offset_z+4];

    // Defects
    b[N_DECISION_VARIABLES+offset_c] = 
        0.5*h[knot]*(problem.knot_params.dx[knot][0] + problem.knot_params.dx[knot+1][0]) + z[offset_z] - z[offset_z+3];
    b[N_DECISION_VARIABLES+offset_c+1] = 
        0.5*h[knot]*(problem.knot_params.dx[knot][1] + problem.knot_params.dx[knot+1][1]) + z[offset_z+1] - z[offset_z+4];

}

void End_Knot_Eq(double *b, double *z, int knot, int offset_z, int offset_c){

}

void Load_Equalities(double *b, double *z, double *h){
    int i, offset_z, offset_c;

    for (i = 0; i < N_KNOTS; i++){
        offset_z = KNOT_SIZE*i;
        Get_Knot_Params(z, offset_z, i);
    }

    offset_z = 0;
    offset_c = 0;
    First_Knot_Eq(b, z, h, 0, offset_z, offset_c);
    offset_c += FIRST_KNOT_CONSTRAINTS;

    for (i = 1; i < N_KNOTS - 1; i++){
        offset_z = i*KNOT_SIZE;
        Middle_Knot_Eq(b, z, h, i, offset_z, offset_c);
        offset_c += MIDDLE_KNOT_CONSTRAINTS;
    }

    offset_z = KNOT_SIZE*(N_KNOTS-1);
    End_Knot_Eq(b, z, N_KNOTS-1, offset_z, offset_c);
}

void First_Knot_Grad(double *b, double *z, int offset){
    b[offset+2] = (mu/(z[offset+2] - ub)) - (mu/(lb - z[offset+2]));
}

void Middle_Knot_Grad(double *b, double *z, int offset){
    b[offset+2] = (mu/(z[offset+2] - ub)) - (mu/(lb - z[offset+2]));
}

void End_Knot_Grad(double *b, double *z, int offset){
    b[offset] = 2.0*(xd[0] - z[offset]);
    b[offset+1] = 2.0*(xd[1] - z[offset+1]);
    b[offset+2] = (mu/(z[offset+2] - ub)) - (mu/(lb - z[offset+2]));
}

void Load_Gradient(double *b, double *z){
    int i, offset;

    offset = 0;
    First_Knot_Grad(b, z, offset);

    for (i = 1; i < N_KNOTS - 1; i++){
        offset = KNOT_SIZE*i;
        Middle_Knot_Grad(b, z, offset);
    }

    offset = KNOT_SIZE*(N_KNOTS-1);
    End_Knot_Grad(b, z, offset);
}

void Load_First_Knot_Columns(int offset_z, int offset_c){
    int i = N_DECISION_VARIABLES;
    int j = i + offset_c;
    int k = offset_z;

    // Jacobian
    problem.k1[0] = Search_For_Sparse_Column(j, k);
    problem.k1[1] = Search_For_Sparse_Column(j, k+1);
    problem.k1[2] = Search_For_Sparse_Column(j, k+3);
    problem.k1[3] = Search_For_Sparse_Column(j, k+4);
    problem.k1[4] = Search_For_Sparse_Column(k, j);
    problem.k1[5] = Search_For_Sparse_Column(k+1, j);
    problem.k1[6] = Search_For_Sparse_Column(k+3, j);
    problem.k1[7] = Search_For_Sparse_Column(k+4, j);

    problem.k1[8] = Search_For_Sparse_Column(j+1, k+1);
    problem.k1[9] = Search_For_Sparse_Column(j+1, k+2);
    problem.k1[10] = Search_For_Sparse_Column(j+1, k+4);
    problem.k1[11] = Search_For_Sparse_Column(j+1, k+5);
    problem.k1[12] = Search_For_Sparse_Column(k+1, j+1);
    problem.k1[13] = Search_For_Sparse_Column(k+2, j+1);
    problem.k1[14] = Search_For_Sparse_Column(k+4, j+1);
    problem.k1[15] = Search_For_Sparse_Column(k+5, j+1);

    problem.k1[16] = Search_For_Sparse_Column(j+2, k);
    problem.k1[17] = Search_For_Sparse_Column(k, j+2);

    problem.k1[18] = Search_For_Sparse_Column(j+3, k+1);
    problem.k1[19] = Search_For_Sparse_Column(k+1, j+3);

    problem.k1[20] = Search_For_Sparse_Column(k+2, k+2);
}

void Load_Middle_Knot_Columns(int knot, int offset_z, int offset_c){
    int i = N_DECISION_VARIABLES;
    int j = i + offset_c;
    int k = offset_z;
    int c = 0;

    problem.km[knot-1][c++] = Search_For_Sparse_Column(j, k);
    problem.km[knot-1][c++] = Search_For_Sparse_Column(j, k+1);
    problem.km[knot-1][c++] = Search_For_Sparse_Column(j, k+3);
    problem.km[knot-1][c++] = Search_For_Sparse_Column(j, k+4);
    problem.km[knot-1][c++] = Search_For_Sparse_Column(k, j);
    problem.km[knot-1][c++] = Search_For_Sparse_Column(k+1, j);
    problem.km[knot-1][c++] = Search_For_Sparse_Column(k+3, j);
    problem.km[knot-1][c++] = Search_For_Sparse_Column(k+4, j);

    problem.km[knot-1][c++] = Search_For_Sparse_Column(j+1, k+1);
    problem.km[knot-1][c++] = Search_For_Sparse_Column(j+1, k+2);
    problem.km[knot-1][c++] = Search_For_Sparse_Column(j+1, k+4);
    problem.km[knot-1][c++] = Search_For_Sparse_Column(j+1, k+5);
    problem.km[knot-1][c++] = Search_For_Sparse_Column(k+1, j+1);
    problem.km[knot-1][c++] = Search_For_Sparse_Column(k+2, j+1);
    problem.km[knot-1][c++] = Search_For_Sparse_Column(k+4, j+1);
    problem.km[knot-1][c++] = Search_For_Sparse_Column(k+5, j+1);

    problem.km[knot-1][c++] = Search_For_Sparse_Column(k+2, k+2);
}

void Load_End_Knot_Columns(int offset_z, int offset_c){
    int k = offset_z;
    int c = 0;

    // hessian
    problem.ke[c++] = Search_For_Sparse_Column(k, k);
    problem.ke[c++] = Search_For_Sparse_Column(k+1, k+1);
    problem.ke[c++] = Search_For_Sparse_Column(k+2, k+2);
}
