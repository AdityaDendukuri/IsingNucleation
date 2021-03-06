#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <stdint.h>
#include <stdbool.h>
#include <time.h>
#include <limits.h>
#include <string.h>
#include "headers/simul.h"
#include "headers/transition_path_sampl.h"


void main()
{
    srand(time(0));
    int* LATTICE;
    int n = 300;
    float J = 1;
    float F = 0.0;
    float T_C = 2.26918531421f;
    float T = 1./0.58; //0.8 * T_C;
    float beta=0.;  
    LATTICE  = malloc(n*n*sizeof(int));
    beta = 1./T;   // J/kb
    
    //run_sim(LATTICE, n, J, F, beta); // experiment 1

    //validate_free_energy(LATTICE, n, J, F); // experiment 2

    F = 0.1;
    //run_sim_nucleation(LATTICE, n, J, F, beta); // experiment 3 with field
    //trans_path_sampl(LATTICE, n,  J,  F,  beta);
    make_nucleus(LATTICE, n, J, F, beta); 
    
    free(LATTICE);
    return;
}
