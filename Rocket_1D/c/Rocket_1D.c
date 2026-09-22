#define _USE_MATH_DEFINES
#include <math.h>
#include <stdlib.h>
#include <stdio.h>
#include "Rocket_1D.h"
#include "SQP.h"

#define MAX(X, Y) ((fabs(X) > (Y)) ? (X) : (Y))
#define A_MIN 0.0

const double ic[3] = {0.0, 0.0, 425.0e3};
const double Shooting_Thrust = 6.0e6;
const double xd = 50.0e3;
const double lb[2] = {25.0e3, 0.0};
const double ub = 8.0e6;
const double tf = 120.0;
const double mu = 1.0e-1;
const double C_param = 3000.0;
const double A_param = 10.52;

typedef struct{
    double dx[N_KNOTS][N_STATES];
    double g[N_KNOTS];
    double rho[N_KNOTS];
    double c_D[N_KNOTS];
    double g_part_h[N_KNOTS];
    double rho_part_h[N_KNOTS];
    double c_D_part_s[N_KNOTS];
    double g_part_h2[N_KNOTS];
    double rho_part_h2[N_KNOTS];
    double c_D_part_s2[N_KNOTS];
}params;

typedef struct{
    int k1[FIRST_KNOT_CONTRIBUTION];
    int km[N_KNOTS-2][MIDDLE_KNOT_CONTRIBUTION];
    int ke[END_KNOT_CONTRIBUTION];
    params knot_params;
}knots;

knots problem;
// knot vars :
// 0 -> altitude,
// 1 -> velocity,
// 2 -> mass,
// 3 -> thrust,
// 4 -> mass slack low,
// 5 -> thrust slack low,
// 6 -> thrust slack high
double Get_tf(){
    return tf;
}

void Load_ic(double *z){

    z[0] = ic[0];
    z[1] = ic[1];
    z[2] = ic[2];

}

void Simulate_Rocket(double *z, double *dt){
    double h, v, m;
    int i, offset_z;

    for (i = 0; i < N_KNOTS; i++){
        offset_z = i*KNOT_SIZE;
        z[offset_z+3] = Shooting_Thrust;
        if (i > 0){
            z[offset_z] = h + dt[i-1]*problem.knot_params.dx[i-1][0];
            z[offset_z+1] = v + dt[i-1]*problem.knot_params.dx[i-1][1];
            z[offset_z+2] = m + dt[i-1]*problem.knot_params.dx[i-1][2];
        }
        Get_Knot_Params(z, offset_z, i);
        h = z[offset_z];
        v = z[offset_z+1];
        m = z[offset_z+2];
    }
}


void Load_slack(double *z){
    int i, j;

    for (i = 0; i < N_KNOTS; i++){
        j = KNOT_SIZE*i;
        z[j+4] = z[j+2] - lb[0];
        z[j+5] = z[j+3] - lb[1];
        z[j+6] = ub - z[j+3];
    }
}

void Load_Problem_Data(double **A, double *b, double *z, double *lambda, double *h, bool sparse){
    // Update and load the gradient of the cost function wrt the decision variables
    Load_Gradient(b, z);

    // Evaluate and load the equality constraints
    // This also updates the knot parameters
    Load_Equalities(b, z, h);

    // Evaluate and load the jacobian of the constraints wrt the decision variables
    Load_Jacobian(A, z, h, sparse);

    // Evaluate and load the hessian of the lagrangian wrt the decision variables
    Load_Hessian(A, z, h, lambda, sparse);

}

double First_Knot_Cost(double *z_knot){
    double contribution = -mu*(log(z_knot[4]) + log(z_knot[5]) + log(z_knot[6]));
    return contribution;
}

double Middle_Knot_Cost(double *z_knot){
    double contribution = -mu*(log(z_knot[4]) + log(z_knot[5]) + log(z_knot[6]));
    return contribution;
}

double End_Knot_Cost(double *z_knot){
    double contribution = pow(z_knot[0] - xd, 2) - 
        mu*(log(z_knot[4]) + log(z_knot[5]) + log(z_knot[6]));
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
    double alpha_check;
    double alpha = 1.0;
    const double tau = 0.995;
    int i, offset_z;
    // First check the maximum alpha that will keep our slack variables positive
    for (i = 0; i < N_DECISION_VARIABLES; i++)
        z_new[i] = z[i] + delta_z[i];

    for (i = 0; i < N_KNOTS; i++){
        offset_z = i*KNOT_SIZE;
        if (z_new[offset_z+4] <= 0.0 && delta_z[offset_z+4] < 0.0){
            alpha_check = - z[offset_z+4] / delta_z[offset_z+4];
            if (alpha_check < alpha)
                alpha = alpha_check;
        }
        if (z_new[offset_z+5] <= 0.0 && delta_z[offset_z+5] < 0.0){
            alpha_check = - z[offset_z+5] / delta_z[offset_z+5];
            if (alpha_check < alpha)
                alpha = alpha_check;
        }
        if (z_new[offset_z+6] <= 0.0 && delta_z[offset_z+6] < 0.0){
            alpha_check = - z[offset_z+6] / delta_z[offset_z+6];
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

void First_Knot_Jacobian(double **A, double *z, double *h, int knot, int offset_z, int offset_c, bool sparse){
    int i = N_DECISION_VARIABLES;
    int j = i + offset_c;
    int k = offset_z;
    const double c1 = -h[knot] / 2.0;
    const double c_A = A_param;
    const double c_C = C_param;

    // const double g1 = problem.knot_params.g[knot];
    const double rho1 = problem.knot_params.rho[knot];
    const double c_D1 = problem.knot_params.c_D[knot];
    const double g_part_h1 = problem.knot_params.g_part_h[knot];
    const double rho_part_h1 = problem.knot_params.rho_part_h[knot];
    const double c_D_part_s1 = problem.knot_params.c_D_part_s[knot];
    // const double h1 = z[offset_z];
    const double v1 = z[offset_z+1];
    const double m1 = z[offset_z+2];
    const double T1 = z[offset_z+3];

    // const double g2 = problem.knot_params.g[knot+1];
    const double rho2 = problem.knot_params.rho[knot+1];
    const double c_D2 = problem.knot_params.c_D[knot+1];
    const double g_part_h2 = problem.knot_params.g_part_h[knot+1];
    const double rho_part_h2 = problem.knot_params.rho_part_h[knot+1];
    const double c_D_part_s2 = problem.knot_params.c_D_part_s[knot+1];
    // const double h2 = z[offset_z+7];
    const double v2 = z[offset_z+8];
    const double m2 = z[offset_z+9];
    const double T2 = z[offset_z+10];

    if (sparse){
        A[j][problem.k1[0]] = -1.0;
        A[j][problem.k1[1]] = c1;
        A[j][problem.k1[2]] = 1.0;
        A[j][problem.k1[3]] = c1;
        A[k][problem.k1[4]] = A[j][problem.k1[0]];
        A[k+1][problem.k1[5]] = A[j][problem.k1[1]];
        A[k+7][problem.k1[6]] = A[j][problem.k1[2]];
        A[k+8][problem.k1[7]] = A[j][problem.k1[3]];

        A[j+1][problem.k1[8]] = MAX(h[knot]*(c_A*v1*v1*c_D1*rho_part_h1 + 2.0*m1*g_part_h1)/(4*m1),A_MIN);
        A[j+1][problem.k1[9]] = MAX((c_A*h[knot]*v1*(v1*c_D_part_s1 + 2.0*c_D1)*rho1/4.0 - m1)/m1,A_MIN);
        A[j+1][problem.k1[10]] = MAX(h[knot]*(-c_A*v1*v1*c_D1*rho1 + 2.0*T1)/(4*m1*m1),A_MIN);
        A[j+1][problem.k1[11]] = c1/m1;
        A[j+1][problem.k1[12]] = MAX(h[knot]*(c_A*v2*v2*c_D2*rho_part_h2 + 2.0*m2*g_part_h2)/(4*m2),A_MIN);
        A[j+1][problem.k1[13]] = MAX((c_A*h[knot]*v2*(v2*c_D_part_s2 + 2.0*c_D2)*rho2/4.0 + m2)/m2,A_MIN);
        A[j+1][problem.k1[14]] = MAX(h[knot]*(-c_A*v2*v2*c_D2*rho2 + 2.0*T2)/(4*m2*m2),A_MIN);
        A[j+1][problem.k1[15]] = c1/m2;
        A[k][problem.k1[16]] = A[j+1][problem.k1[8]];
        A[k+1][problem.k1[17]] = A[j+1][problem.k1[9]];
        A[k+2][problem.k1[18]] = A[j+1][problem.k1[10]];
        A[k+3][problem.k1[19]] = A[j+1][problem.k1[11]];
        A[k+7][problem.k1[20]] = A[j+1][problem.k1[12]];
        A[k+8][problem.k1[21]] = A[j+1][problem.k1[13]];
        A[k+9][problem.k1[22]] = A[j+1][problem.k1[14]];
        A[k+10][problem.k1[23]] = A[j+1][problem.k1[15]]; 

        A[j+2][problem.k1[24]] = -1.0;
        A[j+2][problem.k1[25]] = -c1 / c_C;
        A[j+2][problem.k1[26]] = 1.0;
        A[j+2][problem.k1[27]] = -c1 / c_C;
        A[k+2][problem.k1[28]] = A[j+2][problem.k1[24]];
        A[k+3][problem.k1[29]] = A[j+2][problem.k1[25]];
        A[k+9][problem.k1[30]] = A[j+2][problem.k1[26]];
        A[k+10][problem.k1[31]] = A[j+2][problem.k1[27]];

        A[j+3][problem.k1[32]] = -1.0;
        A[j+3][problem.k1[33]] = 1.0;
        A[k+2][problem.k1[34]] = A[j+3][problem.k1[32]];
        A[k+4][problem.k1[35]] = A[j+3][problem.k1[33]];

        A[j+4][problem.k1[36]] = -1.0;
        A[j+4][problem.k1[37]] = 1.0;
        A[k+3][problem.k1[38]] = A[j+4][problem.k1[36]];
        A[k+5][problem.k1[39]] = A[j+4][problem.k1[37]];

        A[j+5][problem.k1[40]] = 1.0;
        A[j+5][problem.k1[41]] = 1.0;
        A[k+3][problem.k1[42]] = A[j+5][problem.k1[40]];
        A[k+6][problem.k1[43]] = A[j+5][problem.k1[41]];

        A[j+6][problem.k1[44]] = -1.0;
        A[k][problem.k1[45]] = A[j+6][problem.k1[44]];

        A[j+7][problem.k1[46]] = -1.0;
        A[k+1][problem.k1[47]] = A[j+7][problem.k1[46]];

        A[j+8][problem.k1[48]] = -1.0;
        A[k+2][problem.k1[49]] = A[j+8][problem.k1[48]];
    }
    else{
        A[j][k] = -1.0;
        A[j][k+1] = c1;
        A[j][k+7] = 1.0;
        A[j][k+8] = c1;
        A[k][j] = A[j][k];
        A[k+1][j] = A[j][k+1];
        A[k+7][j] = A[j][k+7];
        A[k+8][j] = A[j][k+8];

        A[j+1][k] = MAX(h[knot]*(c_A*v1*v1*c_D1*rho_part_h1 + 2.0*m1*g_part_h1)/(4*m1),A_MIN);
        A[j+1][k+1] = MAX((c_A*h[knot]*v1*(v1*c_D_part_s1 + 2.0*c_D1)*rho1/4.0 - m1)/m1,A_MIN);
        A[j+1][k+2] = MAX(h[knot]*(-c_A*v1*v1*c_D1*rho1 + 2.0*T1)/(4*m1*m1),A_MIN);
        A[j+1][k+3] = MAX(c1/m1,A_MIN);
        A[j+1][k+7] = MAX(h[knot]*(c_A*v2*v2*c_D2*rho_part_h2 + 2.0*m2*g_part_h2)/(4*m2),A_MIN);
        A[j+1][k+8] = MAX((c_A*h[knot]*v2*(v2*c_D_part_s2 + 2.0*c_D2)*rho2/4.0 + m2)/m2,A_MIN);
        A[j+1][k+9] = MAX(h[knot]*(-c_A*v2*v2*c_D2*rho2 + 2.0*T2)/(4*m2*m2),A_MIN);
        A[j+1][k+10] = MAX(c1/m2,A_MIN);
        A[k][j+1] = A[j+1][k];
        A[k+1][j+1] = A[j+1][k+1];
        A[k+2][j+1] = A[j+1][k+2];
        A[k+3][j+1] = A[j+1][k+3];
        A[k+7][j+1] = A[j+1][k+7];
        A[k+8][j+1] = A[j+1][k+8];
        A[k+9][j+1] = A[j+1][k+9];
        A[k+10][j+1] = A[j+1][k+10]; 

        A[j+2][k+2] = -1.0;
        A[j+2][k+3] = -c1 / c_C;
        A[j+2][k+9] = 1.0;
        A[j+2][k+10] = -c1 / c_C;
        A[k+2][j+2] = A[j+2][k+2];
        A[k+3][j+2] = A[j+2][k+3];
        A[k+9][j+2] = A[j+2][k+9];
        A[k+10][j+2] = A[j+2][k+10];

        A[j+3][k+2] = -1.0;
        A[j+3][k+4] = 1.0;
        A[k+2][j+3] = A[j+3][k+2];
        A[k+4][j+3] = A[j+3][k+4];

        A[j+4][k+3] = -1.0;
        A[j+4][k+5] = 1.0;
        A[k+3][j+4] = A[j+4][k+3];
        A[k+5][j+4] = A[j+4][k+5];

        A[j+5][k+3] = 1.0;
        A[j+5][k+6] = 1.0;
        A[k+3][j+5] = A[j+5][k+3];
        A[k+6][j+5] = A[j+5][k+6];

        A[j+6][k] = -1.0;
        A[k][j+6] = A[j+6][k];

        A[j+7][k+1] = -1.0;
        A[k+1][j+7] = A[j+7][k+1];

        A[j+8][k+2] = -1.0;
        A[k+2][j+8] = A[j+8][k+2];
    }
}

void Middle_Knot_Jacobian(double **A, double *z, double *h, int knot, int offset_z, int offset_c, bool sparse){
    int i = N_DECISION_VARIABLES;
    int j = i + offset_c;
    int k = offset_z;
    const double c1 = -h[knot] / 2.0;
    const double c_A = A_param;
    const double c_C = C_param;

    // const double g1 = problem.knot_params.g[knot];
    const double rho1 = problem.knot_params.rho[knot];
    const double c_D1 = problem.knot_params.c_D[knot];
    const double g_part_h1 = problem.knot_params.g_part_h[knot];
    const double rho_part_h1 = problem.knot_params.rho_part_h[knot];
    const double c_D_part_s1 = problem.knot_params.c_D_part_s[knot];
    // const double h1 = z[offset_z];
    const double v1 = z[offset_z+1];
    const double m1 = z[offset_z+2];
    const double T1 = z[offset_z+3];

    // const double g2 = problem.knot_params.g[knot+1];
    const double rho2 = problem.knot_params.rho[knot+1];
    const double c_D2 = problem.knot_params.c_D[knot+1];
    const double g_part_h2 = problem.knot_params.g_part_h[knot+1];
    const double rho_part_h2 = problem.knot_params.rho_part_h[knot+1];
    const double c_D_part_s2 = problem.knot_params.c_D_part_s[knot+1];
    // const double h2 = z[offset_z+7];
    const double v2 = z[offset_z+8];
    const double m2 = z[offset_z+9];
    const double T2 = z[offset_z+10];

    if (sparse){
        A[j][problem.km[knot-1][0]] = -1.0;
        A[j][problem.km[knot-1][1]] = c1;
        A[j][problem.km[knot-1][2]] = 1.0;
        A[j][problem.km[knot-1][3]] = c1;
        A[k][problem.km[knot-1][4]] = A[j][problem.km[knot-1][0]];
        A[k+1][problem.km[knot-1][5]] = A[j][problem.km[knot-1][1]];
        A[k+7][problem.km[knot-1][6]] = A[j][problem.km[knot-1][2]];
        A[k+8][problem.km[knot-1][7]] = A[j][problem.km[knot-1][3]];

        A[j+1][problem.km[knot-1][8]] = MAX(h[knot]*(c_A*v1*v1*c_D1*rho_part_h1 + 2.0*m1*g_part_h1)/(4*m1),A_MIN);
        A[j+1][problem.km[knot-1][9]] = MAX((c_A*h[knot]*v1*(v1*c_D_part_s1 + 2.0*c_D1)*rho1/4.0 - m1)/m1,A_MIN);
        A[j+1][problem.km[knot-1][10]] = MAX(h[knot]*(-c_A*v1*v1*c_D1*rho1 + 2.0*T1)/(4*m1*m1),A_MIN);
        A[j+1][problem.km[knot-1][11]] = c1/m1;
        A[j+1][problem.km[knot-1][12]] = MAX(h[knot]*(c_A*v2*v2*c_D2*rho_part_h2 + 2.0*m2*g_part_h2)/(4*m2),A_MIN);
        A[j+1][problem.km[knot-1][13]] = MAX((c_A*h[knot]*v2*(v2*c_D_part_s2 + 2.0*c_D2)*rho2/4.0 + m2)/m2,A_MIN);
        A[j+1][problem.km[knot-1][14]] = MAX(h[knot]*(-c_A*v2*v2*c_D2*rho2 + 2.0*T2)/(4*m2*m2),A_MIN);
        A[j+1][problem.km[knot-1][15]] = c1/m2;
        A[k][problem.km[knot-1][16]] = A[j+1][problem.km[knot-1][8]];
        A[k+1][problem.km[knot-1][17]] = A[j+1][problem.km[knot-1][9]];
        A[k+2][problem.km[knot-1][18]] = A[j+1][problem.km[knot-1][10]];
        A[k+3][problem.km[knot-1][19]] = A[j+1][problem.km[knot-1][11]];
        A[k+7][problem.km[knot-1][20]] = A[j+1][problem.km[knot-1][12]];
        A[k+8][problem.km[knot-1][21]] = A[j+1][problem.km[knot-1][13]];
        A[k+9][problem.km[knot-1][22]] = A[j+1][problem.km[knot-1][14]];
        A[k+10][problem.km[knot-1][23]] = A[j+1][problem.km[knot-1][15]]; 

        A[j+2][problem.km[knot-1][24]] = -1.0;
        A[j+2][problem.km[knot-1][25]] = -c1 / c_C;
        A[j+2][problem.km[knot-1][26]] = 1.0;
        A[j+2][problem.km[knot-1][27]] = -c1 / c_C;
        A[k+2][problem.km[knot-1][28]] = A[j+2][problem.km[knot-1][24]];
        A[k+3][problem.km[knot-1][29]] = A[j+2][problem.km[knot-1][25]];
        A[k+9][problem.km[knot-1][30]] = A[j+2][problem.km[knot-1][26]];
        A[k+10][problem.km[knot-1][31]] = A[j+2][problem.km[knot-1][27]];

        A[j+3][problem.km[knot-1][32]] = -1.0;
        A[j+3][problem.km[knot-1][33]] = 1.0;
        A[k+2][problem.km[knot-1][34]] = A[j+3][problem.km[knot-1][32]];
        A[k+4][problem.km[knot-1][35]] = A[j+3][problem.km[knot-1][33]];

        A[j+4][problem.km[knot-1][36]] = -1.0;
        A[j+4][problem.km[knot-1][37]] = 1.0;
        A[k+3][problem.km[knot-1][38]] = A[j+4][problem.km[knot-1][36]];
        A[k+5][problem.km[knot-1][39]] = A[j+4][problem.km[knot-1][37]];

        A[j+5][problem.km[knot-1][40]] = 1.0;
        A[j+5][problem.km[knot-1][41]] = 1.0;
        A[k+3][problem.km[knot-1][42]] = A[j+5][problem.km[knot-1][40]];
        A[k+6][problem.km[knot-1][43]] = A[j+5][problem.km[knot-1][41]];
    }
    else{
        A[j][k] = -1.0;
        A[j][k+1] = c1;
        A[j][k+7] = 1.0;
        A[j][k+8] = c1;
        A[k][j] = A[j][k];
        A[k+1][j] = A[j][k+1];
        A[k+7][j] = A[j][k+7];
        A[k+8][j] = A[j][k+8];

        A[j+1][k] = MAX(h[knot]*(c_A*v1*v1*c_D1*rho_part_h1 + 2.0*m1*g_part_h1)/(4*m1),A_MIN);
        A[j+1][k+1] = MAX((c_A*h[knot]*v1*(v1*c_D_part_s1 + 2.0*c_D1)*rho1/4.0 - m1)/m1,A_MIN);
        A[j+1][k+2] = MAX(h[knot]*(-c_A*v1*v1*c_D1*rho1 + 2.0*T1)/(4*m1*m1),A_MIN);
        A[j+1][k+3] = MAX(c1/m1,A_MIN);
        A[j+1][k+7] = MAX(h[knot]*(c_A*v2*v2*c_D2*rho_part_h2 + 2.0*m2*g_part_h2)/(4*m2),A_MIN);
        A[j+1][k+8] = MAX((c_A*h[knot]*v2*(v2*c_D_part_s2 + 2.0*c_D2)*rho2/4.0 + m2)/m2,A_MIN);
        A[j+1][k+9] = MAX(h[knot]*(-c_A*v2*v2*c_D2*rho2 + 2.0*T2)/(4*m2*m2),A_MIN);
        A[j+1][k+10] = MAX(c1/m2,A_MIN);
        A[k][j+1] = A[j+1][k];
        A[k+1][j+1] = A[j+1][k+1];
        A[k+2][j+1] = A[j+1][k+2];
        A[k+3][j+1] = A[j+1][k+3];
        A[k+7][j+1] = A[j+1][k+7];
        A[k+8][j+1] = A[j+1][k+8];
        A[k+9][j+1] = A[j+1][k+9];
        A[k+10][j+1] = A[j+1][k+10]; 

        A[j+2][k+2] = -1.0;
        A[j+2][k+3] = -c1 / c_C;
        A[j+2][k+9] = 1.0;
        A[j+2][k+10] = -c1 / c_C;
        A[k+2][j+2] = A[j+2][k+2];
        A[k+3][j+2] = A[j+2][k+3];
        A[k+9][j+2] = A[j+2][k+9];
        A[k+10][j+2] = A[j+2][k+10];

        A[j+3][k+2] = -1.0;
        A[j+3][k+4] = 1.0;
        A[k+2][j+3] = A[j+3][k+2];
        A[k+4][j+3] = A[j+3][k+4];

        A[j+4][k+3] = -1.0;
        A[j+4][k+5] = 1.0;
        A[k+3][j+4] = A[j+4][k+3];
        A[k+5][j+4] = A[j+4][k+5];

        A[j+5][k+3] = 1.0;
        A[j+5][k+6] = 1.0;
        A[k+3][j+5] = A[j+5][k+3];
        A[k+6][j+5] = A[j+5][k+6];
    }
}

void End_Knot_Jacobian(double **A, int knot, int offset_z, int offset_c, bool sparse){
    int i = N_DECISION_VARIABLES;
    int j = i + offset_c;
    int k = offset_z;

    if (sparse){
        A[j][problem.ke[0]] = -1.0;
        A[j][problem.ke[1]] = 1.0;
        A[k+2][problem.ke[2]] = A[j][problem.ke[0]];
        A[k+4][problem.ke[3]] = A[j][problem.ke[1]];

        A[j+1][problem.ke[4]] = -1.0;
        A[j+1][problem.ke[5]] = 1.0;
        A[k+3][problem.ke[6]] = A[j+1][problem.ke[4]];
        A[k+5][problem.ke[7]] = A[j+1][problem.ke[5]];

        A[j+2][problem.ke[8]] = 1.0;
        A[j+2][problem.ke[9]] = 1.0;
        A[k+3][problem.ke[10]] = A[j+2][problem.ke[8]];
        A[k+6][problem.ke[11]] = A[j+2][problem.ke[9]];
    }
    else{
        A[j][k+2] = -1.0;
        A[j][k+4] = 1.0;
        A[k+2][j] = A[j][k+2];
        A[k+4][j] = A[j][k+4];

        A[j+1][k+3] = -1.0;
        A[j+1][k+5] = 1.0;
        A[k+3][j+1] = A[j+1][k+3];
        A[k+5][j+1] = A[j+1][k+5];

        A[j+2][k+3] = 1.0;
        A[j+2][k+6] = 1.0;
        A[k+3][j+2] = A[j+2][k+3];
        A[k+6][j+2] = A[j+2][k+6];
    }
}

void Load_Jacobian(double **A, double *z, double *h, bool sparse){
    int i, offset_z, offset_c;
    offset_z = 0;
    offset_c = 0;

    First_Knot_Jacobian(A, z, h, 0, offset_z, offset_c, sparse);
    offset_c += FIRST_KNOT_CONSTRAINTS;

    for (i = 1; i < N_KNOTS - 1; i++){
        offset_z = i*KNOT_SIZE;
        Middle_Knot_Jacobian(A, z, h, i, offset_z, offset_c, sparse);
        offset_c += MIDDLE_KNOT_CONSTRAINTS;
    }

    offset_z = KNOT_SIZE*(N_KNOTS-1);
    End_Knot_Jacobian(A, N_KNOTS-1, offset_z, offset_c, sparse);
}

void First_Knot_Hessian(double **A, double *z, double *h, double *lambda,
                        int knot, int offset_z, int offset_c, bool sparse){
    int i = offset_z;
    const double c_A = A_param;
    // const double c_C = C_param;

    // const double g1 = problem.knot_params.g[knot];
    const double rho1 = problem.knot_params.rho[knot];
    const double c_D1 = problem.knot_params.c_D[knot];
    // const double g_part_h = problem.knot_params.g_part_h[knot];
    const double rho_part_h = problem.knot_params.rho_part_h[knot];
    const double c_D_part_s = problem.knot_params.c_D_part_s[knot];
    const double g_part_h2 = problem.knot_params.g_part_h2[knot];
    const double rho_part_h2 = problem.knot_params.rho_part_h2[knot];
    const double c_D_part_s2 = problem.knot_params.c_D_part_s2[knot];
    // const double h1 = z[offset_z];
    const double v1 = z[offset_z+1];
    const double m1 = z[offset_z+2];
    const double T1 = z[offset_z+3];
    const double s_m_lb = z[offset_z+4];
    const double s_T_lb = z[offset_z+5];
    const double s_T_ub = z[offset_z+6];
    const double lambda2 = lambda[offset_c+1];

    if (sparse){
        A[i][problem.k1[50]] = MAX((h[knot]*lambda2*(c_A*v1*v1*c_D1*rho_part_h2 + 2.0*m1*g_part_h2))/(4*m1),A_MIN);
        A[i][problem.k1[51]] = MAX((c_A*h[knot]*lambda2*v1*(v1*c_D_part_s + 2.0*c_D1)*rho_part_h)/(4*m1),A_MIN);
        A[i][problem.k1[52]] = MAX((-c_A*h[knot]*lambda2*v1*v1*c_D1*rho_part_h)/(4*m1*m1),A_MIN);
        
        A[i+1][problem.k1[53]] = A[i][problem.k1[51]];
        A[i+1][problem.k1[54]] = MAX(c_A*h[knot]*lambda2*(v1*v1*c_D_part_s2 + 4.0*v1*c_D_part_s + 2.0*c_D1)*rho1/(4.0*m1),A_MIN);
        A[i+1][problem.k1[55]] = MAX(-c_A*h[knot]*lambda2*v1*(v1*c_D_part_s + 2.0*c_D1)*rho1/(4.0*m1*m1),A_MIN);

        A[i+2][problem.k1[56]] = A[i][problem.k1[52]];
        A[i+2][problem.k1[57]] = A[i+1][problem.k1[55]];
        A[i+2][problem.k1[58]] = MAX(h[knot]*lambda2*(c_A*v1*v1*c_D1*rho1 - 2.0*T1)/(2.0*m1*m1*m1),A_MIN);
        A[i+2][problem.k1[59]] = MAX(h[knot]*lambda2/(2.0*m1*m1),A_MIN);

        A[i+3][problem.k1[60]] = A[i+2][problem.k1[59]];
        
        A[i+4][problem.k1[61]] = mu / pow(s_m_lb,2);

        A[i+5][problem.k1[62]] = mu / pow(s_T_lb,2);

        A[i+6][problem.k1[63]] = mu / pow(s_T_ub,2);
    }
    else{
        A[i][i] = MAX((h[knot]*lambda2*(c_A*v1*v1*c_D1*rho_part_h2 + 2.0*m1*g_part_h2))/(4*m1),A_MIN);
        A[i][i+1] = MAX((c_A*h[knot]*lambda2*v1*(v1*c_D_part_s + 2.0*c_D1)*rho_part_h)/(4*m1),A_MIN);
        A[i][i+2] = MAX((-c_A*h[knot]*lambda2*v1*v1*c_D1*rho_part_h)/(4*m1*m1),A_MIN);
        
        A[i+1][i] = A[i][i+1];
        A[i+1][i+1] = MAX(c_A*h[knot]*lambda2*(v1*v1*c_D_part_s2 + 4.0*v1*c_D_part_s + 2.0*c_D1)*rho1/(4.0*m1),A_MIN);
        A[i+1][i+2] = MAX(-c_A*h[knot]*lambda2*v1*(v1*c_D_part_s + 2.0*c_D1)*rho1/(4.0*m1*m1),A_MIN);

        A[i+2][i] = A[i][i+2];
        A[i+2][i+1] = A[i+1][i+2];
        A[i+2][i+2] = MAX(h[knot]*lambda2*(c_A*v1*v1*c_D1*rho1 - 2.0*T1)/(2.0*m1*m1*m1),A_MIN);
        A[i+2][i+3] = MAX(h[knot]*lambda2/(2.0*m1*m1),A_MIN);

        A[i+3][i+2] = A[i+2][i+3];
        
        A[i+4][i+4] = mu / pow(s_m_lb,2);

        A[i+5][i+5] = mu / pow(s_T_lb,2);

        A[i+6][i+6] = mu / pow(s_T_ub,2);
    }
}

void Middle_Knot_Hessian(double **A, double *z, double *h, double *lambda,
                         int knot, int offset_z, int offset_c, int offset_c_last, bool sparse){
    int i = offset_z;
    const double c_A = A_param;
    // const double c_C = C_param;

    // const double g1 = problem.knot_params.g[knot];
    const double rho1 = problem.knot_params.rho[knot];
    const double c_D1 = problem.knot_params.c_D[knot];
    // const double g_part_h = problem.knot_params.g_part_h[knot];
    const double rho_part_h = problem.knot_params.rho_part_h[knot];
    const double c_D_part_s = problem.knot_params.c_D_part_s[knot];
    const double g_part_h2 = problem.knot_params.g_part_h2[knot];
    const double rho_part_h2 = problem.knot_params.rho_part_h2[knot];
    const double c_D_part_s2 = problem.knot_params.c_D_part_s2[knot];
    // const double h1 = z[offset_z];
    const double v1 = z[offset_z+1];
    const double m1 = z[offset_z+2];
    const double T1 = z[offset_z+3];
    const double s_m_lb = z[offset_z+4];
    const double s_T_lb = z[offset_z+5];
    const double s_T_ub = z[offset_z+6];
    const double lambda2 = lambda[offset_c+1] + lambda[offset_c_last+1];

    if (sparse){
        A[i][problem.km[knot-1][44]] = MAX((h[knot]*lambda2*(c_A*v1*v1*c_D1*rho_part_h2 + 2.0*m1*g_part_h2))/(4*m1),A_MIN);
        A[i][problem.km[knot-1][45]] = MAX((c_A*h[knot]*lambda2*v1*(v1*c_D_part_s + 2.0*c_D1)*rho_part_h)/(4*m1),A_MIN);
        A[i][problem.km[knot-1][46]] = MAX((-c_A*h[knot]*lambda2*v1*v1*c_D1*rho_part_h)/(4*m1*m1),A_MIN);
        
        A[i+1][problem.km[knot-1][47]] = A[i][problem.km[knot-1][45]];
        A[i+1][problem.km[knot-1][48]] = MAX(c_A*h[knot]*lambda2*(v1*v1*c_D_part_s2 + 4.0*v1*c_D_part_s + 2.0*c_D1)*rho1/(4.0*m1),A_MIN);
        A[i+1][problem.km[knot-1][49]] = MAX(-c_A*h[knot]*lambda2*v1*(v1*c_D_part_s + 2.0*c_D1)*rho1/(4.0*m1*m1),A_MIN);

        A[i+2][problem.km[knot-1][50]] = A[i][problem.km[knot-1][46]];
        A[i+2][problem.km[knot-1][51]] = A[i+1][problem.km[knot-1][49]];
        A[i+2][problem.km[knot-1][52]] = MAX(h[knot]*lambda2*(c_A*v1*v1*c_D1*rho1 - 2.0*T1)/(2.0*m1*m1*m1),A_MIN);
        A[i+2][problem.km[knot-1][53]] = MAX(h[knot]*lambda2/(2.0*m1*m1),A_MIN);

        A[i+3][problem.km[knot-1][54]] = A[i+2][problem.km[knot-1][53]];
        
        A[i+4][problem.km[knot-1][55]] = mu / pow(s_m_lb,2);

        A[i+5][problem.km[knot-1][56]] = mu / pow(s_T_lb,2);

        A[i+6][problem.km[knot-1][57]] = mu / pow(s_T_ub,2);
    }
    else{
        A[i][i] = MAX((h[knot]*lambda2*(c_A*v1*v1*c_D1*rho_part_h2 + 2.0*m1*g_part_h2))/(4*m1),A_MIN);
        A[i][i+1] = MAX((c_A*h[knot]*lambda2*v1*(v1*c_D_part_s + 2.0*c_D1)*rho_part_h)/(4*m1),A_MIN);
        A[i][i+2] = MAX((-c_A*h[knot]*lambda2*v1*v1*c_D1*rho_part_h)/(4*m1*m1),A_MIN);
        
        A[i+1][i] = A[i][i+1];
        A[i+1][i+1] = MAX(c_A*h[knot]*lambda2*(v1*v1*c_D_part_s2 + 4.0*v1*c_D_part_s + 2.0*c_D1)*rho1/(4.0*m1),A_MIN);
        A[i+1][i+2] = MAX(-c_A*h[knot]*lambda2*v1*(v1*c_D_part_s + 2.0*c_D1)*rho1/(4.0*m1*m1),A_MIN);

        A[i+2][i] = A[i][i+2];
        A[i+2][i+1] = A[i+1][i+2];
        A[i+2][i+2] = MAX(h[knot]*lambda2*(c_A*v1*v1*c_D1*rho1 - 2.0*T1)/(2.0*m1*m1*m1),A_MIN);
        A[i+2][i+3] = MAX(h[knot]*lambda2/(2.0*m1*m1),A_MIN);

        A[i+3][i+2] = A[i+2][i+3];
        
        A[i+4][i+4] = mu / pow(s_m_lb,2);

        A[i+5][i+5] = mu / pow(s_T_lb,2);

        A[i+6][i+6] = mu / pow(s_T_ub,2);
    }
}

void End_Knot_Hessian(double **A, double *z, double *h, double *lambda,
                      int knot, int offset_z, int offset_c, bool sparse){
    int i = offset_z;
    const double c_A = A_param;
    // const double c_C = C_param;

    // const double g1 = problem.knot_params.g[knot];
    const double rho1 = problem.knot_params.rho[knot];
    const double c_D1 = problem.knot_params.c_D[knot];
    // const double g_part_h = problem.knot_params.g_part_h[knot];
    const double rho_part_h = problem.knot_params.rho_part_h[knot];
    const double c_D_part_s = problem.knot_params.c_D_part_s[knot];
    const double g_part_h2 = problem.knot_params.g_part_h2[knot];
    const double rho_part_h2 = problem.knot_params.rho_part_h2[knot];
    const double c_D_part_s2 = problem.knot_params.c_D_part_s2[knot];
    // const double h1 = z[offset_z];
    const double v1 = z[offset_z+1];
    const double m1 = z[offset_z+2];
    const double T1 = z[offset_z+3];
    const double s_m_lb = z[offset_z+4];
    const double s_T_lb = z[offset_z+5];
    const double s_T_ub = z[offset_z+6];
    const double lambda2 = lambda[offset_c+1];

    if (sparse){
        A[i][problem.ke[12]] = MAX(2.0 + (h[knot]*lambda2*(c_A*v1*v1*c_D1*rho_part_h2 + 2.0*m1*g_part_h2))/(4*m1),A_MIN);
        A[i][problem.ke[13]] = MAX((c_A*h[knot]*lambda2*v1*(v1*c_D_part_s + 2.0*c_D1)*rho_part_h)/(4*m1),A_MIN);
        A[i][problem.ke[14]] = MAX((-c_A*h[knot]*lambda2*v1*v1*c_D1*rho_part_h)/(4*m1*m1),A_MIN);
        
        A[i+1][problem.ke[15]] = A[i][problem.ke[13]];
        A[i+1][problem.ke[16]] = MAX(c_A*h[knot]*lambda2*(v1*v1*c_D_part_s2 + 4.0*v1*c_D_part_s + 2.0*c_D1)*rho1/(4.0*m1),A_MIN);
        A[i+1][problem.ke[17]] = MAX(-c_A*h[knot]*lambda2*v1*(v1*c_D_part_s + 2.0*c_D1)*rho1/(4.0*m1*m1),A_MIN);

        A[i+2][problem.ke[18]] = A[i][problem.ke[14]];
        A[i+2][problem.ke[19]] = A[i+1][problem.ke[17]];
        A[i+2][problem.ke[20]] = MAX(h[knot]*lambda2*(c_A*v1*v1*c_D1*rho1 - 2.0*T1)/(2.0*m1*m1*m1),A_MIN);
        A[i+2][problem.ke[21]] = MAX(h[knot]*lambda2/(2.0*m1*m1),A_MIN);

        A[i+3][problem.ke[22]] = A[i+2][problem.ke[21]];
        
        A[i+4][problem.ke[23]] = mu / pow(s_m_lb,2);

        A[i+5][problem.ke[24]] = mu / pow(s_T_lb,2);

        A[i+6][problem.ke[25]] = mu / pow(s_T_ub,2);
    }
    else{
        A[i][i] = MAX(2.0 + (h[knot]*lambda2*(c_A*v1*v1*c_D1*rho_part_h2 + 2.0*m1*g_part_h2))/(4*m1),A_MIN);
        A[i][i+1] = MAX((c_A*h[knot]*lambda2*v1*(v1*c_D_part_s + 2.0*c_D1)*rho_part_h)/(4*m1),A_MIN);
        A[i][i+2] = MAX((-c_A*h[knot]*lambda2*v1*v1*c_D1*rho_part_h)/(4*m1*m1),A_MIN);
        
        A[i+1][i] = A[i][i+1];
        A[i+1][i+1] = MAX(c_A*h[knot]*lambda2*(v1*v1*c_D_part_s2 + 4.0*v1*c_D_part_s + 2.0*c_D1)*rho1/(4.0*m1),A_MIN);
        A[i+1][i+2] = MAX(-c_A*h[knot]*lambda2*v1*(v1*c_D_part_s + 2.0*c_D1)*rho1/(4.0*m1*m1),A_MIN);

        A[i+2][i] = A[i][i+2];
        A[i+2][i+1] = A[i+1][i+2];
        A[i+2][i+2] = MAX(h[knot]*lambda2*(c_A*v1*v1*c_D1*rho1 - 2.0*T1)/(2.0*m1*m1*m1),A_MIN);
        A[i+2][i+3] = MAX(h[knot]*lambda2/(2.0*m1*m1),A_MIN);

        A[i+3][i+2] = A[i+2][i+3];
        
        A[i+4][i+4] = mu / pow(s_m_lb,2);

        A[i+5][i+5] = mu / pow(s_T_lb,2);

        A[i+6][i+6] = mu / pow(s_T_ub,2);
    }

}

void Load_Hessian(double **A, double *z, double *h, double *lambda, bool sparse){
    int i, offset_z, offset_c, offset_c_last;
    offset_z = 0;
    offset_c = 0;

    First_Knot_Hessian(A, z, h, lambda, 0, offset_z, offset_c, sparse);
    offset_c_last = offset_c;
    offset_c += FIRST_KNOT_CONSTRAINTS;

    for (i = 1; i < N_KNOTS - 1; i++){
        offset_z = i*KNOT_SIZE;
        Middle_Knot_Hessian(A, z, h, lambda, i, offset_z, offset_c, offset_c_last, sparse);
        offset_c_last = offset_c;
        offset_c += MIDDLE_KNOT_CONSTRAINTS;
    }

    offset_z = KNOT_SIZE*(N_KNOTS-1);
    End_Knot_Hessian(A, z, h, lambda, N_KNOTS-1, offset_z, offset_c, sparse);
}

void Get_Knot_Params(double *z, int offset, int knot){
    double c1, c2, c3, c4, c5, c6, c7;
    double Temperature, Pressure, Speed_of_Sound, Mach, Drag;
    const double gravity0 = 9.8065;
    const double Re = 6371.0e3;

    // Extract knot variables
    double h = z[offset];
    double v = z[offset+1];
    double m = z[offset+2];
    double T = z[offset+3];

    // Gravity calculation
    
    problem.knot_params.g[knot] = gravity0*pow(Re/(Re + h), 2);
    problem.knot_params.g_part_h[knot] = -2.0*Re*gravity0/pow(Re + h, 3);
    problem.knot_params.g_part_h2[knot] = (6.0*Re*gravity0)/pow(Re + h, 4);

    // Air density calculation
    if (h < 11000.0){
        c1 = 15.04;
        c2 = 0.00649;
        c3 = 101.29;
        c4 = 273.1;
        c5 = 288.08;
        c6 = 5.256;
        c7 = 0.2869;

        Temperature = c1 - c2*h;
        Pressure =c3*pow((Temperature+c4)/c5, c6);

        problem.knot_params.rho_part_h[knot] = (c2*c3*pow((c1 - c2*h + c4)/c5, c6)*
                             (1 - c6))/(c7*pow(c1-c2*h+c4, 2));

        problem.knot_params.rho_part_h2[knot] = (pow(c2, 2)*c3*pow((c1 - c2*h + c4)/c5, c6)*
                             (c6 - 1)*(c6 - 2))/(c7*pow(c1-c2*h+c4, 3));
    }
    else if (h < 25000.0){
        c1 = -56.46;
        c2 = 22.65;
        c3 = 1.73;
        c4 = 0.000157;
        c5 = 0.2869;
        c6 = 273.1;
        Temperature = c1;
        Pressure = c2*exp(c3 - c4*h);

        problem.knot_params.rho_part_h[knot] = (-c2*c4*exp(c3 - c4*h))/(c5*(c1+c6));
        problem.knot_params.rho_part_h2[knot] = (c2*pow(c4,2)*exp(c3 - c4*h))/(c5*(c1+c6));
    }
    else{
        c1 = -131.21;
        c2 = 0.00299;
        c3 = 2.488;
        c4 = 273.1;
        c5 = 216.6;
        c6 = -11.388;
        c7 = 0.2869;
        Temperature = c1 + c2*h;
        Pressure = c3*(pow((Temperature+c4)/c5, c6));

        problem.knot_params.rho_part_h[knot] = (c2*c3*pow((c1 + c2*h + c4)/c5, c6)*
                                (c6 - 1))/(c7*pow(c1+c2*h+c4, 2));
        problem.knot_params.rho_part_h2[knot] = (pow(c2, 2)*c3*pow((c1 + c2*h + c4)/c5, c6)*
                                (c6 - 1)*(c6 - 2))/(c7*pow(c1+c2*h+c4, 3));
    }
    
    problem.knot_params.rho[knot] = Pressure/(0.2869*(Temperature+273.1));

    Speed_of_Sound = sqrt(1.4*287*(Temperature+273.1));
    Mach = v / Speed_of_Sound;

    if (Mach < 0.8){
        problem.knot_params.c_D[knot] = 0.22;
        problem.knot_params.c_D_part_s[knot] = 0.0;
        problem.knot_params.c_D_part_s2[knot] = 0.0;
    }
    else if (Mach < 1.2){
        c1 = 0.22;
        c2 = 0.48;
        c3 = 0.8;

        problem.knot_params.c_D[knot] = c1 + c2*pow(sin(M_PI*(Mach - c3)/c3), 2);
        problem.knot_params.c_D_part_s[knot] = 
                (M_PI*c2*sin(2*M_PI*v/(Speed_of_Sound*c3)))/(Speed_of_Sound*c3);
        problem.knot_params.c_D_part_s2[knot] = 
                (2*pow(M_PI, 2)*c2*cos((2*M_PI*v)/(Speed_of_Sound*c3)))/(pow(Speed_of_Sound,2)*pow(c3,2));
    }
    else{
        c1 = 0.25;
        c2 = 0.54;
        c3 = 1.2;

        problem.knot_params.c_D[knot] = c1 + c2/pow(Mach,c3);
        problem.knot_params.c_D_part_s[knot] = 
                (-c2*c3*pow(v/Speed_of_Sound, -c3))/v;
        problem.knot_params.c_D_part_s2[knot] = 
                (c2*c3*pow(v/Speed_of_Sound, -c3)*(c3 + 1))/pow(v,2);
    }

    Drag = 0.5 * A_param * problem.knot_params.c_D[knot] * problem.knot_params.rho[knot] * pow(v, 2);

    problem.knot_params.dx[knot][0] = v;
    problem.knot_params.dx[knot][1] = (T - Drag)/m - problem.knot_params.g[knot];
    problem.knot_params.dx[knot][2] = -T / C_param;
}

void First_Knot_Eq(double *b, double *z, double *h, int knot, int offset_z, int offset_c){
    // h1 = z[offset_z]
    // v1 = z[offset_z+1]
    // m1 = z[offset_z+2]
    // T1 = z[offset_z+3]
    // s_m_lb = z[offset_z+4]
    // s_T_lb = z[offset_z+5]
    // s_T_ub = z[offset_z+6]
    // h2 = z[offset_z+7]
    // v2 = z[offset_z+8]
    // m2 = z[offset_z+9]
    // T2 = z[offset_z+10]

    // Defects
    b[N_DECISION_VARIABLES+offset_c] = 
        0.5*h[knot]*(problem.knot_params.dx[knot][0] + problem.knot_params.dx[knot+1][0]) + z[offset_z] - z[offset_z+7];
    b[N_DECISION_VARIABLES+offset_c+1] = 
        0.5*h[knot]*(problem.knot_params.dx[knot][1] + problem.knot_params.dx[knot+1][1]) + z[offset_z+1] - z[offset_z+8];
    b[N_DECISION_VARIABLES+offset_c+2] = 
        0.5*h[knot]*(problem.knot_params.dx[knot][2] + problem.knot_params.dx[knot+1][2]) + z[offset_z+2] - z[offset_z+9];
    // Lower bounds
    b[N_DECISION_VARIABLES+offset_c+3] = z[offset_z+2] - lb[0] - z[offset_z+4];
    b[N_DECISION_VARIABLES+offset_c+4] = z[offset_z+3] - lb[1] - z[offset_z+5];
    // Upper bounds
    b[N_DECISION_VARIABLES+offset_c+5] = ub - z[offset_z+3] - z[offset_z+6];
    // Initial conditions
    b[N_DECISION_VARIABLES+offset_c+6] = z[offset_z] - ic[0];
    b[N_DECISION_VARIABLES+offset_c+7] = z[offset_z+1] - ic[1];
    b[N_DECISION_VARIABLES+offset_c+8] = z[offset_z+2] - ic[2];
}

void Middle_Knot_Eq(double *b, double *z, double *h, int knot, int offset_z, int offset_c){
    // h1 = z[offset_z]
    // v1 = z[offset_z+1]
    // m1 = z[offset_z+2]
    // T1 = z[offset_z+3]
    // s_m_lb = z[offset_z+4]
    // s_T_lb = z[offset_z+5]
    // s_T_ub = z[offset_z+6]
    // h2 = z[offset_z+7]
    // v2 = z[offset_z+8]
    // m2 = z[offset_z+9]
    // T2 = z[offset_z+10]

    // Defects
    b[N_DECISION_VARIABLES+offset_c] = 
        0.5*h[knot]*(problem.knot_params.dx[knot][0] + problem.knot_params.dx[knot+1][0]) + z[offset_z] - z[offset_z+7];
    b[N_DECISION_VARIABLES+offset_c+1] = 
        0.5*h[knot]*(problem.knot_params.dx[knot][1] + problem.knot_params.dx[knot+1][1]) + z[offset_z+1] - z[offset_z+8];
    b[N_DECISION_VARIABLES+offset_c+2] = 
        0.5*h[knot]*(problem.knot_params.dx[knot][2] + problem.knot_params.dx[knot+1][2]) + z[offset_z+2] - z[offset_z+9];
    // Lower bounds
    b[N_DECISION_VARIABLES+offset_c+3] = z[offset_z+2] - lb[0] - z[offset_z+4];
    b[N_DECISION_VARIABLES+offset_c+4] = z[offset_z+3] - lb[1] - z[offset_z+5];
    // Upper bounds
    b[N_DECISION_VARIABLES+offset_c+5] = ub - z[offset_z+3] - z[offset_z+6];

}

void End_Knot_Eq(double *b, double *z, int knot, int offset_z, int offset_c){
    // h1 = z[offset_z]
    // v1 = z[offset_z+1]
    // m1 = z[offset_z+2]
    // T1 = z[offset_z+3]
    // s_m_lb = z[offset_z+4]
    // s_T_lb = z[offset_z+5]
    // s_T_ub = z[offset_z+6]

    // Lower bounds
    b[N_DECISION_VARIABLES+offset_c] = z[offset_z+2] - lb[0] - z[offset_z+4];
    b[N_DECISION_VARIABLES+offset_c+1] = z[offset_z+3] - lb[1] - z[offset_z+5];
    // Upper bounds
    b[N_DECISION_VARIABLES+offset_c+2] = ub - z[offset_z+3] - z[offset_z+6];
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
    b[offset+4] = mu / z[offset+4];
    b[offset+5] = mu / z[offset+5];
    b[offset+6] = mu / z[offset+6];
}

void Middle_Knot_Grad(double *b, double *z, int offset){
    b[offset+4] = mu / z[offset+4];
    b[offset+5] = mu / z[offset+5];
    b[offset+6] = mu / z[offset+6];
}

void End_Knot_Grad(double *b, double *z, int offset){
    b[offset] = 2.0*(xd - z[offset]);
    b[offset+4] = mu / z[offset+4];
    b[offset+5] = mu / z[offset+5];
    b[offset+6] = mu / z[offset+6];
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

void Load_First_Knot_Columns(int offset_z, int offset_c, int **Checklist_S, bool sparse){
    int i = N_DECISION_VARIABLES;
    int j = i + offset_c;
    int k = offset_z;
    int c = 0;
    const int rows[FIRST_KNOT_CONTRIBUTION] = {
        //jacobian
        j, j, j, j, k, k+1, k+7, k+8,
        j+1, j+1, j+1, j+1, j+1, j+1, j+1, j+1, k, k+1, k+2, k+3, k+7, k+8, k+9, k+10,
        j+2, j+2, j+2, j+2, k+2, k+3, k+9, k+10,
        j+3, j+3, k+2, k+4,
        j+4, j+4, k+3, k+5,
        j+5, j+5, k+3, k+6,
        j+6, k,
        j+7, k+1,
        j+8, k+2,
        // hessian
        k, k, k,
        k+1, k+1, k+1,
        k+2, k+2, k+2, k+2,
        k+3,
        k+4,
        k+5,
        k+6
    };
    const int cols[FIRST_KNOT_CONTRIBUTION] = {
        //jacobian
        k, k+1, k+7, k+8, j, j, j, j,
        k, k+1, k+2, k+3, k+7, k+8, k+9, k+10, j+1, j+1, j+1, j+1, j+1, j+1, j+1, j+1,
        k+2, k+3, k+9, k+10, j+2, j+2, j+2, j+2,
        k+2, k+4, j+3, j+3,
        k+3, k+5, j+4, j+4,
        k+3, k+6, j+5, j+5,
        k, j+6,
        k+1, j+7,
        k+2, j+8,
        // hessian
        k, k+1, k+2,
        k, k+1, k+2,
        k, k+1, k+2, k+3,
        k+2,
        k+4,
        k+5,
        k+6
    };
    for (c = 0; c < FIRST_KNOT_CONTRIBUTION; c++)
        if (sparse)
            Checklist_S[rows[c]][cols[c]] = 1;
        else
            problem.k1[c] = Search_For_Sparse_Column(rows[c],cols[c]);

}

void Load_Middle_Knot_Columns(int knot, int offset_z, int offset_c, int **Checklist_S, bool sparse){
    int i = N_DECISION_VARIABLES;
    int j = i + offset_c;
    int k = offset_z;
    int c = 0;
    const int rows[MIDDLE_KNOT_CONTRIBUTION] = {
        //jacobian
        j, j, j, j, k, k+1, k+7, k+8,
        j+1, j+1, j+1, j+1, j+1, j+1, j+1, j+1, k, k+1, k+2, k+3, k+7, k+8, k+9, k+10,
        j+2, j+2, j+2, j+2, k+2, k+3, k+9, k+10,
        j+3, j+3, k+2, k+4,
        j+4, j+4, k+3, k+5,
        j+5, j+5, k+3, k+6,
        // hessian
        k, k, k,
        k+1, k+1, k+1,
        k+2, k+2, k+2, k+2,
        k+3,
        k+4,
        k+5,
        k+6
    };
    const int cols[MIDDLE_KNOT_CONTRIBUTION] = {
        //jacobian
        k, k+1, k+7, k+8, j, j, j, j,
        k, k+1, k+2, k+3, k+7, k+8, k+9, k+10, j+1, j+1, j+1, j+1, j+1, j+1, j+1, j+1,
        k+2, k+3, k+9, k+10, j+2, j+2, j+2, j+2,
        k+2, k+4, j+3, j+3,
        k+3, k+5, j+4, j+4,
        k+3, k+6, j+5, j+5,
        // hessian
        k, k+1, k+2,
        k, k+1, k+2,
        k, k+1, k+2, k+3,
        k+2,
        k+4,
        k+5,
        k+6
    };
    for (c = 0; c < MIDDLE_KNOT_CONTRIBUTION; c++)
        if (sparse)
            Checklist_S[rows[c]][cols[c]] = 1;
        else
            problem.km[knot-1][c] = Search_For_Sparse_Column(rows[c],cols[c]);
   
}

void Load_End_Knot_Columns(int offset_z, int offset_c, int **Checklist_S, bool sparse){
    int i = N_DECISION_VARIABLES;
    int j = i + offset_c;
    int k = offset_z;
    int c = 0;

    const int rows[END_KNOT_CONTRIBUTION] = {
        //jacobian
        j, j, k+2, k+4,
        j+1, j+1, k+3, k+5,
        j+2, j+2, k+3, k+6,
        // hessian
        k, k, k,
        k+1, k+1, k+1,
        k+2, k+2, k+2, k+2,
        k+3,
        k+4,
        k+5,
        k+6
    };
    const int cols[END_KNOT_CONTRIBUTION] = {
        // jacobian
        k+2, k+4, j, j,
        k+3, k+5, j+1, j+1,
        k+3, k+6, j+2, j+2,
        // hessian
        k, k+1, k+2,
        k, k+1, k+2,
        k, k+1, k+2, k+3,
        k+2,
        k+4,
        k+5,
        k+6
    };

    for (c = 0; c < END_KNOT_CONTRIBUTION; c++)
        if (sparse)
            Checklist_S[rows[c]][cols[c]] = 1;
        else
            problem.ke[c] = Search_For_Sparse_Column(rows[c],cols[c]);

}
