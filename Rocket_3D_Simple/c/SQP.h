#ifndef SQP_H
#define SQP_H

#include <time.h>
#include <stdbool.h>

typedef struct{
    double pf;
    double df;
    int iterations;
    bool converged;
    struct timespec time_load, time_elim, time_subs, time_extract, time_alpha, time_update;
}sqp_results;

void Initial_Problem_Data();

void Initial_Solve();

void Cleanup();

sqp_results Run_SQP();

double Primary_Feasability_Check();

double Dual_Feasability_Check();

void Reset_A_S();

int Search_For_Sparse_Column(int row, int col);

void Get_IC_Cols(bool sparse);

void Parse_Checklist_S();

void Fill_A_S_Zeros();

void Convert_id_U_S();

double *Get_z();

#endif