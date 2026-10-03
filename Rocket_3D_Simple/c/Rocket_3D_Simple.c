#define _USE_MATH_DEFINES
#include <math.h>
#include <stdlib.h>
#include <stdio.h>
#include "Rocket_3D_Simple.h"
#include "SQP.h"

#define MAX(X, Y) ((fabs(X) > (Y)) ? (X) : (Y))
#define A_MIN 0.0

const double initial_fpa = M_PI/4;
const double initial_hdg = 0.0;
const double Shooting_c_L_z = -0.1;
const double ic[8] = {0.0, 0.0, -100.0, 50.0, 0.0, 0.0, 0.0, 0.0};
const double xd[3] = {800.0, 200.0, -300.0};
const double c_L_lb = -1.5;
const double c_L_ub = 1.5;
const double p_d_ub = 0.0;
const double tf = 10.0;
const double mu = 1.0e-2;
const double mu_q = 1.0e1;
const double A_param = 0.25;

typedef struct{
    double t[N_KNOTS];
    double h[N_KNOTS-1];
    double dx[N_KNOTS][N_STATES];
    double T[N_KNOTS];
    double m[N_KNOTS];
    double g[N_KNOTS];
    double rho[N_KNOTS];
    double c_D[N_KNOTS];
    double rho_part_d[N_KNOTS];
    double c_D_part_s[N_KNOTS];
    double rho_part_d2[N_KNOTS];
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
// 0 -> North position,
// 1 -> East position,
// 2 -> Down position,
// 3 -> Speed,
// 4 -> NED to velocity quaternion scalar
// 5 -> NED to velocity quaternion imaginary x
// 6 -> NED to velocity quaternion imaginary y
// 7 -> NED to velocity quaternion imaginary z
// 8 -> wy
// 9 -> wz
// 10 -> c_Ly
// 11 -> c_Lz

double Get_tf(){
    return tf;
}

double *Get_t(){
    return &problem.knot_params.t[0];
}

void Set_t(){
    double dt = tf / (double)(N_KNOTS - 1);
    problem.knot_params.t[0] = 0.0;
    for (int i = 1; i < N_KNOTS; i++){
        problem.knot_params.t[i] =problem.knot_params.t[i-1] + dt;
        problem.knot_params.h[i-1] = dt; 
    }
}

void Load_ic(double *z){
    const double gam_2 = initial_fpa / 2.0;
    const double chi_2 = initial_hdg / 2.0;
    z[0] = ic[0];
    z[1] = ic[1];
    z[2] = ic[2];
    z[3] = ic[3];
    z[4] = cos(gam_2)*cos(chi_2);
    z[5] = -sin(gam_2)*sin(chi_2);
    z[6] = sin(gam_2)*cos(chi_2);
    z[7] = cos(gam_2)*sin(chi_2);

}

void Simulate_Rocket(double *z){
    double p_N, p_E, p_D, s, q0, q1, q2, q3, wy, wz;
    double delta_q[4], theta, w_norm, dt;
    int i, offset_z;

    for (i = 0; i < N_KNOTS; i++){
        offset_z = i*KNOT_SIZE;
        z[offset_z+10] = 0.0;
        z[offset_z+11] = Shooting_c_L_z;
        if (i > 0){
            dt = problem.knot_params.h[i-1];
            z[offset_z] = p_N + dt*problem.knot_params.dx[i-1][0];
            z[offset_z+1] = p_E + dt*problem.knot_params.dx[i-1][1];
            z[offset_z+2] = p_D + dt*problem.knot_params.dx[i-1][2];
            z[offset_z+3] = s + dt*problem.knot_params.dx[i-1][3];
            w_norm = sqrt(pow(wy,2) + pow(wz,2));
            theta = w_norm*dt;
            delta_q[0] = cos(theta/2.0);
            delta_q[2] = (wy/w_norm)*sin(theta/2.0);
            delta_q[3] = (wz/w_norm)*sin(theta/2.0);
            z[offset_z+4] = q0*delta_q[0] - q2*delta_q[2] - q3*delta_q[3];
            z[offset_z+5] = q1*delta_q[0] - q3*delta_q[2] + q2*delta_q[3];
            z[offset_z+6] = q2*delta_q[0] + q0*delta_q[2] - q1*delta_q[3];
            z[offset_z+7] = q3*delta_q[0] + q1*delta_q[2] + q0*delta_q[3];
        }
        // Make sure to set w in this func
        Get_Knot_Params(z, offset_z, i);
        p_N = z[offset_z];
        p_E = z[offset_z+1];
        p_D = z[offset_z+2];
        s   = z[offset_z+3];
        q0  = z[offset_z+4];
        q1  = z[offset_z+5];
        q2  = z[offset_z+6];
        q3  = z[offset_z+7];
        wy  = problem.knot_params.dx[i][8];
        wz  = problem.knot_params.dx[i][9];
    }
}

void Load_Problem_Data(double **A, double *b, double *z, double *lambda, bool sparse){
    // Update and load the gradient of the cost function wrt the decision variables
    Load_Gradient(b, z);

    // Evaluate and load the equality constraints
    // This also updates the knot parameters
    Load_Equalities(b, z);

    // Evaluate and load the jacobian of the constraints wrt the decision variables
    Load_Jacobian(A, z, sparse);

    // Evaluate and load the hessian of the lagrangian wrt the decision variables
    Load_Hessian(A, z, lambda, sparse);

}

double First_Knot_Cost(double *z_knot){
    double contribution;
    // double contribution = -mu*(log(z_knot[2] - lb[0]) +
    //                            log(z_knot[3] - lb[1]) +
    //                            log(ub - z_knot[3]));
    return contribution;
}

double Middle_Knot_Cost(double *z_knot){
    double contribution;
    // double contribution = -mu*(log(z_knot[2] - lb[0]) +
    //                            log(z_knot[3] - lb[1]) +
    //                            log(ub - z_knot[3]));
    return contribution;
}

double End_Knot_Cost(double *z_knot){
    double contribution;
    // double contribution = pow(z_knot[0] - xd, 2) - 
    //                     mu*(log(z_knot[2] - lb[0]) +
    //                         log(z_knot[3] - lb[1]) +
    //                         log(ub - z_knot[3]));
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
    // First check the maximum alpha that will keep our slack variables positive
    for (i = 0; i < N_DECISION_VARIABLES; i++)
        z_new[i] = z[i] + delta_z[i];

    // for (i = 0; i < N_KNOTS; i++){
    //     offset_z = i*KNOT_SIZE;
    //     slack_new[0] = z_new[offset_z+2] - lb[0];
    //     slack_new[1] = z_new[offset_z+3] - lb[1];
    //     slack_new[2] = ub - z_new[offset_z+3];
    //     slack_cur[0] = z[offset_z+2] - lb[0];
    //     slack_cur[1] = z[offset_z+3] - lb[1];
    //     slack_cur[2] = ub - z[offset_z+3];
    //     for (j = 0; j < N_LBS + N_UBS; j++){
    //         if (slack_new[j] <= 0.0) {
    //             alpha_check = -slack_cur[j] / (slack_new[j] - slack_cur[j]);
    //             if (alpha_check < alpha)
    //                 alpha = alpha_check;
    //         }
    //     }

    // }

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

void First_Knot_Jacobian(double **A, double *z, int knot, int offset_z, int offset_c, bool sparse){
    int i = N_DECISION_VARIABLES;
    int j = i + offset_c;
    int k = offset_z;
    const double c1 = -problem.knot_params.h[knot] / 2.0;
    const double c_A = A_param;

    // const double g1 = problem.knot_params.g[knot];
    const double rho1 = problem.knot_params.rho[knot];
    const double c_D1 = problem.knot_params.c_D[knot];
    const double rho_part_d1 = problem.knot_params.rho_part_d[knot];
    const double c_D_part_s1 = problem.knot_params.c_D_part_s[knot];
    // const double h1 = z[offset_z];
    const double v1 = z[offset_z+1];
    const double m1 = z[offset_z+2];
    const double T1 = z[offset_z+3];

    // const double g2 = problem.knot_params.g[knot+1];
    const double rho2 = problem.knot_params.rho[knot+1];
    const double c_D2 = problem.knot_params.c_D[knot+1];
    const double rho_part_d2 = problem.knot_params.rho_part_d[knot+1];
    const double c_D_part_s2 = problem.knot_params.c_D_part_s[knot+1];
    // const double h2 = z[offset_z+4];
    const double v2 = z[offset_z+5];
    const double m2 = z[offset_z+6];
    const double T2 = z[offset_z+7];

    if (sparse){
        
    }
    else{
        A[j][k] = -1.0;
        A[j][k+1] = c1;
        A[j][k+4] = 1.0;
        A[j][k+5] = c1;
        A[k][j] = A[j][k];
        A[k+1][j] = A[j][k+1];
        A[k+4][j] = A[j][k+4];
        A[k+5][j] = A[j][k+5];

    }
}

void Middle_Knot_Jacobian(double **A, double *z, int knot, int offset_z, int offset_c, bool sparse){
    int i = N_DECISION_VARIABLES;
    int j = i + offset_c;
    int k = offset_z;
    const double c1 = -problem.knot_params.h[knot] / 2.0;
    const double c_A = A_param;

    // const double g1 = problem.knot_params.g[knot];
    const double rho1 = problem.knot_params.rho[knot];
    const double c_D1 = problem.knot_params.c_D[knot];
    const double rho_part_d1 = problem.knot_params.rho_part_d[knot];
    const double c_D_part_s1 = problem.knot_params.c_D_part_s[knot];
    // const double h1 = z[offset_z];
    const double v1 = z[offset_z+1];
    const double m1 = z[offset_z+2];
    const double T1 = z[offset_z+3];

    // const double g2 = problem.knot_params.g[knot+1];
    const double rho2 = problem.knot_params.rho[knot+1];
    const double c_D2 = problem.knot_params.c_D[knot+1];
    const double rho_part_d2 = problem.knot_params.rho_part_d[knot+1];
    const double c_D_part_s2 = problem.knot_params.c_D_part_s[knot+1];
    // const double h2 = z[offset_z+4];
    const double v2 = z[offset_z+5];
    const double m2 = z[offset_z+6];
    const double T2 = z[offset_z+7];

    if (sparse){

    }
    else{
        A[j][k] = -1.0;
        A[j][k+1] = c1;
        A[j][k+4] = 1.0;
        A[j][k+5] = c1;
        A[k][j] = A[j][k];
        A[k+1][j] = A[j][k+1];
        A[k+4][j] = A[j][k+4];
        A[k+5][j] = A[j][k+5];

    }
}

void End_Knot_Jacobian(double **A, int knot, int offset_z, int offset_c, bool sparse){
    
}

void Load_Jacobian(double **A, double *z, bool sparse){
    int i, offset_z, offset_c;
    offset_z = 0;
    offset_c = 0;

    First_Knot_Jacobian(A, z, 0, offset_z, offset_c, sparse);
    offset_c += FIRST_KNOT_CONSTRAINTS;

    for (i = 1; i < N_KNOTS - 1; i++){
        offset_z = i*KNOT_SIZE;
        Middle_Knot_Jacobian(A, z, i, offset_z, offset_c, sparse);
        offset_c += MIDDLE_KNOT_CONSTRAINTS;
    }

    offset_z = KNOT_SIZE*(N_KNOTS-1);
    End_Knot_Jacobian(A, N_KNOTS-1, offset_z, offset_c, sparse);
}

void First_Knot_Hessian(double **A, double *z, double *lambda,
                        int knot, int offset_z, int offset_c, bool sparse){
    int i = offset_z;
    const double c_A = A_param;
    // const double c_C = C_param;

    // const double g1 = problem.knot_params.g[knot];
    const double rho1 = problem.knot_params.rho[knot];
    const double c_D1 = problem.knot_params.c_D[knot];
    const double rho_part_d = problem.knot_params.rho_part_d[knot];
    const double c_D_part_s = problem.knot_params.c_D_part_s[knot];
    const double rho_part_d2 = problem.knot_params.rho_part_d2[knot];
    const double c_D_part_s2 = problem.knot_params.c_D_part_s2[knot];
    // const double h1 = z[offset_z];
    const double v1 = z[offset_z+1];
    const double m1 = z[offset_z+2];
    const double T1 = z[offset_z+3];

    const double lambda2 = lambda[offset_c+1];

    if (sparse){
       
    }
    else{

    }
}

void Middle_Knot_Hessian(double **A, double *z, double *lambda,
                         int knot, int offset_z, int offset_c, int offset_c_last, bool sparse){
    int i = offset_z;
    const double c_A = A_param;
    // const double c_C = C_param;

    // const double g1 = problem.knot_params.g[knot];
    const double rho1 = problem.knot_params.rho[knot];
    const double c_D1 = problem.knot_params.c_D[knot];
    const double rho_part_d = problem.knot_params.rho_part_d[knot];
    const double c_D_part_s = problem.knot_params.c_D_part_s[knot];
    const double rho_part_d2 = problem.knot_params.rho_part_d2[knot];
    const double c_D_part_s2 = problem.knot_params.c_D_part_s2[knot];
    // const double h1 = z[offset_z];
    const double v1 = z[offset_z+1];
    const double m1 = z[offset_z+2];
    const double T1 = z[offset_z+3];

    const double lambda2 = lambda[offset_c+1] + lambda[offset_c_last+1];

    if (sparse){ 

    }
    else{

    }
}

void End_Knot_Hessian(double **A, double *z, double *lambda,
                      int knot, int offset_z, int offset_c, bool sparse){
    int i = offset_z;
    const double c_A = A_param;
    // const double c_C = C_param;

    // const double g1 = problem.knot_params.g[knot];
    const double rho1 = problem.knot_params.rho[knot];
    const double c_D1 = problem.knot_params.c_D[knot];
    const double rho_part_d = problem.knot_params.rho_part_d[knot];
    const double c_D_part_s = problem.knot_params.c_D_part_s[knot];
    const double rho_part_d2 = problem.knot_params.rho_part_d2[knot];
    const double c_D_part_s2 = problem.knot_params.c_D_part_s2[knot];
    // const double h1 = z[offset_z];
    const double v1 = z[offset_z+1];
    const double m1 = z[offset_z+2];
    const double T1 = z[offset_z+3];

    const double lambda2 = lambda[offset_c+1];

    if (sparse){

    }
    else{
        
    }

}

void Load_Hessian(double **A, double *z, double *lambda, bool sparse){
    int i, offset_z, offset_c, offset_c_last;
    offset_z = 0;
    offset_c = 0;

    First_Knot_Hessian(A, z, lambda, 0, offset_z, offset_c, sparse);
    offset_c_last = offset_c;
    offset_c += FIRST_KNOT_CONSTRAINTS;

    for (i = 1; i < N_KNOTS - 1; i++){
        offset_z = i*KNOT_SIZE;
        Middle_Knot_Hessian(A, z, lambda, i, offset_z, offset_c, offset_c_last, sparse);
        offset_c_last = offset_c;
        offset_c += MIDDLE_KNOT_CONSTRAINTS;
    }

    offset_z = KNOT_SIZE*(N_KNOTS-1);
    End_Knot_Hessian(A, z, lambda, N_KNOTS-1, offset_z, offset_c, sparse);
}

void Get_Knot_Params(double *z, int offset, int knot){
    double c1, c2, c3, c4, c5, c6, c7;
    double ay, az;
    double Temperature, Pressure, Speed_of_Sound, Mach, Drag;
    const double gravity0 = 9.8065;
    const double Re = 6371.0e3;

    // Extract knot variables
    double p_N = z[offset];
    double p_E = z[offset+1];
    double p_D = z[offset+2];
    double s = z[offset+3];
    double q0 = z[offset+4];
    double q1 = z[offset+5];
    double q2 = z[offset+6];
    double q3 = z[offset+7];
    double wy = z[offset+8];
    double wz = z[offset+9];
    double c_Ly = z[offset+10];
    double c_Lz = z[offset+11];

    // Gravity calculation
    
    problem.knot_params.g[knot] = gravity0*pow(Re/(Re - p_D), 2);

    // Air density calculation
    if (p_D > -11000.0){
        c1 = 15.04;
        c2 = 0.00649;
        c3 = 101.29;
        c4 = 273.1;
        c5 = 288.08;
        c6 = 5.256;
        c7 = 0.2869;

        Temperature = c1 + c2*p_D;
        Pressure =c3*pow((Temperature+c4)/c5, c6);

        problem.knot_params.rho_part_d[knot] = (c2*c3*pow((c1 + c2*p_D + c4)/c5, c6)*
                             (c6 - 1.0))/(c7*pow(c1+c2*p_D+c4, 2));

        problem.knot_params.rho_part_d2[knot] = (pow(c2, 2)*c3*pow((c1 + c2*p_D + c4)/c5, c6)*
                             (pow(c6,2) - 3.0*c6 + 2.0))/(c7*pow(c1+c2*p_D+c4, 3));
    }
    else if (p_D >= -25000.0){
        c1 = -56.46;
        c2 = 22.65;
        c3 = 1.73;
        c4 = 0.000157;
        c5 = 0.2869;
        c6 = 273.1;
        Temperature = c1;
        Pressure = c2*exp(c3 + c4*p_D);

        problem.knot_params.rho_part_d[knot] = (c2*c4*exp(c3 + c4*p_D))/(c5*(c1+c6));
        problem.knot_params.rho_part_d2[knot] = (c2*pow(c4,2)*exp(c3 + c4*p_D))/(c5*(c1+c6));
    }
    else{
        c1 = -131.21;
        c2 = 0.00299;
        c3 = 2.488;
        c4 = 273.1;
        c5 = 216.6;
        c6 = -11.388;
        c7 = 0.2869;
        Temperature = c1 - c2*p_D;
        Pressure = c3*(pow((Temperature+c4)/c5, c6));

        problem.knot_params.rho_part_d[knot] = (c2*c3*pow((c1 - c2*p_D + c4)/c5, c6)*
                                (1.0 - c6))/(c7*pow(c1-c2*p_D+c4, 2));
        problem.knot_params.rho_part_d2[knot] = (pow(c2, 2)*c3*pow((c1 - c2*p_D + c4)/c5, c6)*
                                (pow(c6, 2) - 3.0*c6 + 2.0))/(c7*pow(c1-c2*p_D+c4, 3));
    }
    
    problem.knot_params.rho[knot] = Pressure/(0.2869*(Temperature+273.1));

    Speed_of_Sound = sqrt(1.4*287*(Temperature+273.1));
    Mach = s / Speed_of_Sound;

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
                (M_PI*c2*sin(2*M_PI*s/(Speed_of_Sound*c3)))/(Speed_of_Sound*c3);
        problem.knot_params.c_D_part_s2[knot] = 
                (2*pow(M_PI, 2)*c2*cos((2*M_PI*s)/(Speed_of_Sound*c3)))/(pow(Speed_of_Sound,2)*pow(c3,2));
    }
    else{
        c1 = 0.25;
        c2 = 0.54;
        c3 = 1.2;

        problem.knot_params.c_D[knot] = c1 + c2/pow(Mach,c3);
        problem.knot_params.c_D_part_s[knot] = 
                (-c2*c3*pow(s/Speed_of_Sound, -c3))/s;
        problem.knot_params.c_D_part_s2[knot] = 
                (c2*c3*pow(s/Speed_of_Sound, -c3)*(c3 + 1))/pow(s,2);
    }

    const double t_table[] = {0.0, 0.2, 0.5, 2.5, 3.0, 3.25, 4.0, 6.0, 8.0, 10.0, 11.0, 12.0, 13.0, 13.5};
    const double T_table[] = {0.0, 300.0, 1000.0, 1000.0, 800.0, 600.0, 550.0, 525.0, 500.0, 450.0, 350.0, 250.0, 100.0, 0.0};
    const double m_table[] = {15, 14.92, 14.52, 11.8533, 11.32, 11.12, 10.57, 9.17, 7.8367, 6.6367, 6.17, 5.8367, 5.7033, 5.7033};
    problem.knot_params.T[knot] = T_table[(sizeof(T_table)/sizeof(T_table[0]))-1];
    problem.knot_params.m[knot] = m_table[(sizeof(m_table)/sizeof(m_table[0]))-1];
    for (int i = 1; i < sizeof(t_table)/sizeof(t_table[0]); i++){
        if (problem.knot_params.t[knot] <= t_table[i]){
            c1 = (problem.knot_params.t[knot] - t_table[i-1])/(t_table[i] - t_table[i-1]);
            c2 = 1.0 - c1;
            problem.knot_params.T[knot] = c1*T_table[i] + c2*T_table[i-1];
            problem.knot_params.m[knot] = c1*m_table[i] + c2*m_table[i-1];
            break;
        } 
    }

    Drag = 0.5 * A_param * problem.knot_params.c_D[knot] * problem.knot_params.rho[knot] * pow(s, 2);
    az = 0.5 * A_param * c_Lz * problem.knot_params.rho[knot] * pow(s,2) / problem.knot_params.m[knot];
    ay = 0.5 * A_param * c_Ly * problem.knot_params.rho[knot] * pow(s,2) / problem.knot_params.m[knot];

    problem.knot_params.dx[knot][0] = s*(-2.0*pow(q2,2) - 2.0*pow(q3,2) + 1.0);
    problem.knot_params.dx[knot][1] = 2.0*s*(q0*q3 + q1*q2);
    problem.knot_params.dx[knot][2] = 2.0*s*(-q0*q2 + q1*q3);

    problem.knot_params.dx[knot][3] = ((problem.knot_params.T[knot] - Drag)/problem.knot_params.m[knot]) + 2.0*problem.knot_params.g[knot]*(-q0*q2 + q1*q3);

    problem.knot_params.dx[knot][4] = 0.5*(-wy*q2 - wz*q3);
    problem.knot_params.dx[knot][5] = 0.5*(wz*q2 - wy*q3);
    problem.knot_params.dx[knot][6] = 0.5*(wy*q0 - wz*q1);
    problem.knot_params.dx[knot][7] = 0.5*(wz*q0 + wy*q1);

    problem.knot_params.dx[knot][8] = (-az - problem.knot_params.g[knot]*(-2.0*pow(q1,2) - 2.0*pow(q2,2) + 1.0)) / s;
    problem.knot_params.dx[knot][9] = (ay + 2.0*problem.knot_params.g[knot]*(q0*q1 + q2*q3)) / s;

}

void First_Knot_Eq(double *b, double *z, int knot, int offset_z, int offset_c){
    // h1 = z[offset_z]
    // v1 = z[offset_z+1]
    // m1 = z[offset_z+2]
    // T1 = z[offset_z+3]
    // h2 = z[offset_z+4]
    // v2 = z[offset_z+5]
    // m2 = z[offset_z+6]
    // T2 = z[offset_z+7]

    // Defects
    // b[N_DECISION_VARIABLES+offset_c] = 
    //     0.5*h[knot]*(problem.knot_params.dx[knot][0] + problem.knot_params.dx[knot+1][0]) + z[offset_z] - z[offset_z+4];
    // b[N_DECISION_VARIABLES+offset_c+1] = 
    //     0.5*h[knot]*(problem.knot_params.dx[knot][1] + problem.knot_params.dx[knot+1][1]) + z[offset_z+1] - z[offset_z+5];
    // b[N_DECISION_VARIABLES+offset_c+2] = 
    //     0.5*h[knot]*(problem.knot_params.dx[knot][2] + problem.knot_params.dx[knot+1][2]) + z[offset_z+2] - z[offset_z+6];
    // // Initial conditions
    // b[N_DECISION_VARIABLES+offset_c+3] = z[offset_z] - ic[0];
    // b[N_DECISION_VARIABLES+offset_c+4] = z[offset_z+1] - ic[1];
    // b[N_DECISION_VARIABLES+offset_c+5] = z[offset_z+2] - ic[2];
}

void Middle_Knot_Eq(double *b, double *z, int knot, int offset_z, int offset_c){
    // h1 = z[offset_z]
    // v1 = z[offset_z+1]
    // m1 = z[offset_z+2]
    // T1 = z[offset_z+3]
    // h2 = z[offset_z+4]
    // v2 = z[offset_z+5]
    // m2 = z[offset_z+6]
    // T2 = z[offset_z+7]

    // Defects
    // b[N_DECISION_VARIABLES+offset_c] = 
    //     0.5*h[knot]*(problem.knot_params.dx[knot][0] + problem.knot_params.dx[knot+1][0]) + z[offset_z] - z[offset_z+4];
    // b[N_DECISION_VARIABLES+offset_c+1] = 
    //     0.5*h[knot]*(problem.knot_params.dx[knot][1] + problem.knot_params.dx[knot+1][1]) + z[offset_z+1] - z[offset_z+5];
    // b[N_DECISION_VARIABLES+offset_c+2] = 
    //     0.5*h[knot]*(problem.knot_params.dx[knot][2] + problem.knot_params.dx[knot+1][2]) + z[offset_z+2] - z[offset_z+6];

}

void End_Knot_Eq(double *b, double *z, int knot, int offset_z, int offset_c){

}

void Load_Equalities(double *b, double *z){
    int i, offset_z, offset_c;

    for (i = 0; i < N_KNOTS; i++){
        offset_z = KNOT_SIZE*i;
        Get_Knot_Params(z, offset_z, i);
    }

    offset_z = 0;
    offset_c = 0;
    First_Knot_Eq(b, z, 0, offset_z, offset_c);
    offset_c += FIRST_KNOT_CONSTRAINTS;

    for (i = 1; i < N_KNOTS - 1; i++){
        offset_z = i*KNOT_SIZE;
        Middle_Knot_Eq(b, z, i, offset_z, offset_c);
        offset_c += MIDDLE_KNOT_CONSTRAINTS;
    }

    offset_z = KNOT_SIZE*(N_KNOTS-1);
    End_Knot_Eq(b, z, N_KNOTS-1, offset_z, offset_c);
}

void First_Knot_Grad(double *b, double *z, int offset){
    // b[offset+2] = -mu/(lb[0] - z[offset+2]);
    // b[offset+3] = -mu/(lb[1] - z[offset+3]) + mu/(z[offset+3] - ub);
}

void Middle_Knot_Grad(double *b, double *z, int offset){
    // b[offset+2] = -mu/(lb[0] - z[offset+2]);
    // b[offset+3] = -mu/(lb[1] - z[offset+3]) + mu/(z[offset+3] - ub);
}

void End_Knot_Grad(double *b, double *z, int offset){
    // b[offset] = 2.0*(xd - z[offset]);
    // b[offset+2] = -mu/(lb[0] - z[offset+2]);
    // b[offset+3] = -mu/(lb[1] - z[offset+3]) + mu/(z[offset+3] - ub);
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
        j, j, j, j, k, k+1, k+4, k+5,
        j+1, j+1, j+1, j+1, j+1, j+1, j+1, j+1, k, k+1, k+2, k+3, k+4, k+5, k+6, k+7,
        j+2, j+2, j+2, j+2, k+2, k+3, k+6, k+7,
        j+3, k,
        j+4, k+1,
        j+5, k+2,
        // hessian
        k, k, k,
        k+1, k+1, k+1,
        k+2, k+2, k+2, k+2,
        k+3, k+3
    };
    const int cols[FIRST_KNOT_CONTRIBUTION] = {
        //jacobian
        k, k+1, k+4, k+5, j, j, j, j,
        k, k+1, k+2, k+3, k+4, k+5, k+6, k+7, j+1, j+1, j+1, j+1, j+1, j+1, j+1, j+1,
        k+2, k+3, k+6, k+7, j+2, j+2, j+2, j+2,
        k, j+3,
        k+1, j+4,
        k+2, j+5,
        // hessian
        k, k+1, k+2,
        k, k+1, k+2,
        k, k+1, k+2, k+3,
        k+2, k+3
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
        j, j, j, j, k, k+1, k+4, k+5,
        j+1, j+1, j+1, j+1, j+1, j+1, j+1, j+1, k, k+1, k+2, k+3, k+4, k+5, k+6, k+7,
        j+2, j+2, j+2, j+2, k+2, k+3, k+6, k+7,
        // hessian
        k, k, k,
        k+1, k+1, k+1,
        k+2, k+2, k+2, k+2,
        k+3, k+3
    };
    const int cols[MIDDLE_KNOT_CONTRIBUTION] = {
        //jacobian
        k, k+1, k+4, k+5, j, j, j, j,
        k, k+1, k+2, k+3, k+4, k+5, k+6, k+7, j+1, j+1, j+1, j+1, j+1, j+1, j+1, j+1,
        k+2, k+3, k+6, k+7, j+2, j+2, j+2, j+2,
        // hessian
        k, k+1, k+2,
        k, k+1, k+2,
        k, k+1, k+2, k+3,
        k+2, k+3
    };
    for (c = 0; c < MIDDLE_KNOT_CONTRIBUTION; c++)
        if (sparse)
            Checklist_S[rows[c]][cols[c]] = 1;
        else
            problem.km[knot-1][c] = Search_For_Sparse_Column(rows[c],cols[c]);
   
}

void Load_End_Knot_Columns(int offset_z, int offset_c, int **Checklist_S, bool sparse){
    int k = offset_z;
    int c = 0;

    const int rows[END_KNOT_CONTRIBUTION] = {
        //jacobian
        // hessian
        k, k, k,
        k+1, k+1, k+1,
        k+2, k+2, k+2, k+2,
        k+3, k+3
    };
    const int cols[END_KNOT_CONTRIBUTION] = {
        // jacobian
        // hessian
        k, k+1, k+2,
        k, k+1, k+2,
        k, k+1, k+2, k+3,
        k+2, k+3
    };

    for (c = 0; c < END_KNOT_CONTRIBUTION; c++)
        if (sparse)
            Checklist_S[rows[c]][cols[c]] = 1;
        else
            problem.ke[c] = Search_For_Sparse_Column(rows[c],cols[c]);

}
