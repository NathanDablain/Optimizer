#ifndef ROCKET_3D_SIMPLE_H
#define ROCKET_3D_SIMPLE_H

#include <stdbool.h>

#define N_STATES 10
#define N_INPUTS 2
#define N_KNOTS 10
#define N_LBS 2
#define N_UBS 3
#define FIRST_KNOT_CONSTRAINTS 18
#define MIDDLE_KNOT_CONSTRAINTS 10
#define END_KNOT_CONSTRAINTS 2
#define FIRST_KNOT_CONTRIBUTION (2*(19) + 12)
#define MIDDLE_KNOT_CONTRIBUTION (2*(16) + 12)
#define END_KNOT_CONTRIBUTION (2*(0) + 12)
#define TOTAL_CONTRIBUTION (FIRST_KNOT_CONTRIBUTION + MIDDLE_KNOT_CONTRIBUTION*(N_KNOTS-2) + END_KNOT_CONTRIBUTION)
#define KNOT_SIZE (N_STATES + N_INPUTS)
#define N_DECISION_VARIABLES (N_KNOTS * KNOT_SIZE)
#define N_CONSTRAINTS (FIRST_KNOT_CONSTRAINTS + MIDDLE_KNOT_CONSTRAINTS*(N_KNOTS-2) + END_KNOT_CONSTRAINTS)
#define SYSTEM_SIZE (N_DECISION_VARIABLES + N_CONSTRAINTS)

double Get_tf();

void Set_t();

void Load_ic(double *z);

void Simulate_Rocket(double *z);

void Load_Problem_Data(double **A, double *b, double *z, double *lambda, bool sparse);

double Select_Alpha(double *z, double *delta_z);

double Get_Cost(double *z);

void Load_Jacobian(double **A, double *z, bool sparse);

void Load_Hessian(double **A, double *z, double *lambda, bool sparse);

void Get_Knot_Params(double *z, int offset, int knot);

void Load_Equalities(double *b, double *z);

void Load_Gradient(double *b, double *z);

void Load_First_Knot_Columns(int offset_z, int offset_c, int **Checklist_S, bool sparse);

void Load_Middle_Knot_Columns(int knot, int offset_z, int offset_c, int **Checklist_S, bool sparse);

void Load_End_Knot_Columns(int offset_z, int offset_c, int **Checklist_S, bool sparse);

double *Get_t();

#endif