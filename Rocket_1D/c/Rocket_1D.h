#ifndef ROCKET_1D_H
#define ROCKET_1D_H

#include <stdbool.h>

#define N_STATES 3
#define N_INPUTS 1
#define N_KNOTS 100
#define N_LBS 2
#define N_UBS 1
#define FIRST_KNOT_CONSTRAINTS 6
#define MIDDLE_KNOT_CONSTRAINTS 3
#define END_KNOT_CONSTRAINTS 0
#define FIRST_KNOT_CONTRIBUTION (2*(19) + 12)
#define MIDDLE_KNOT_CONTRIBUTION (2*(16) + 12)
#define END_KNOT_CONTRIBUTION (2*(0) + 12)
#define TOTAL_CONTRIBUTION (FIRST_KNOT_CONTRIBUTION + MIDDLE_KNOT_CONTRIBUTION*(N_KNOTS-2) + END_KNOT_CONTRIBUTION)
#define KNOT_SIZE (N_STATES + N_INPUTS)
#define N_DECISION_VARIABLES (N_KNOTS * KNOT_SIZE)
#define N_CONSTRAINTS (FIRST_KNOT_CONSTRAINTS + MIDDLE_KNOT_CONSTRAINTS*(N_KNOTS-2) + END_KNOT_CONSTRAINTS)
#define SYSTEM_SIZE (N_DECISION_VARIABLES + N_CONSTRAINTS)

double Get_tf();

void Load_ic(double *z);

void Simulate_Rocket(double *z, double *h);

void Load_Problem_Data(double **A, double *b, double *z, double *lambda, double *h, bool sparse);

double Select_Alpha(double *z, double *delta_z);

double Get_Cost(double *z);

void Load_Jacobian(double **A, double *z, double *h, bool sparse);

void Load_Hessian(double **A, double *z, double *h, double *lambda, bool sparse);

void Get_Knot_Params(double *z, int offset, int knot);

void Load_Equalities(double *b, double *z, double *h);

void Load_Gradient(double *b, double *z);

void Load_First_Knot_Columns(int offset_z, int offset_c, int **Checklist_S, bool sparse);

void Load_Middle_Knot_Columns(int knot, int offset_z, int offset_c, int **Checklist_S, bool sparse);

void Load_End_Knot_Columns(int offset_z, int offset_c, int **Checklist_S, bool sparse);

#endif