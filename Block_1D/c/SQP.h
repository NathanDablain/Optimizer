#ifndef SQP_H
#define SQP_H

#include <stdbool.h>

typedef struct{
    double pf;
    double df;
    int iterations;
    bool converged;
}sqp_results;

void Initial_Problem_Data();

void Initial_Solve();

void Cleanup();

sqp_results Run_SQP();

double Primary_Feasability_Check();

double Dual_Feasability_Check();

void Reset_A_S();

int Search_For_Sparse_Column(int row, int col);

void Get_IC_Cols();

void Parse_Checklist_S();

void Fill_A_S_Zeros();

void Convert_id_U_S();

double *Get_z();

double *Get_t();

#endif