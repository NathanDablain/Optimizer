#include <stdint.h>
#include <stdlib.h>
#include <stdio.h>
#include <math.h>
#include <time.h>
#include <string.h>
#include "LU.h"
#include "Block_1D.h"
#include "SQP.h"

int plot_flag = 1;

struct timespec time_load, time_elim, time_subs, time_extract, time_alpha, time_update;

int main(){
    int i;
    struct timespec time_start, time_end;

    Initial_Problem_Data();

    Initial_Solve();

    timespec_get(&time_start, TIME_UTC);

    // Run problem
    sqp_results results = Run_SQP();

    timespec_get(&time_end, TIME_UTC);

    Cleanup();

    // Log data
    double *z = Get_z();
    double *t = Get_t();

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

    printf("\n\n");
    printf("%.9f seconds elapsed\n", time_elapsed);
    if (results.converged)
        printf("Converged in %d iterations\n", results.iterations);
    else
        printf("Failed to converge in %d iterations\n", results.iterations);
    printf("Primal feasability: %.5e\n", results.pf);
    printf("Dual feasability: %.5e\n", results.df);

    return 1;
}
