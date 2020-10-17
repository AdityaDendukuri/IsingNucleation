#include "../headers/monte_carlo.h"
#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <stdint.h>
#include <stdbool.h>
#include <time.h>
#include <limits.h>
#include <string.h>
#include "../headers/util.h"
#include "../headers/ising_lattice.h"


// single flip monte carlo step (boltzmann)
void mc_step_boltz(int* LATTICE, int n,  float beta, float J, float F) { 
    int x, y = 0;
    int num_spins = n*n;
    for(int i = 0; i < num_spins; i++) { 
        x = (int) rand()%n;
        y = (int) rand()%n;
        // hamiltonian
        float delE =  2 * J * ( LATTICE[ y + n*PBC(x+1, n) ] +
                                LATTICE[ PBC(y+1, n) + n*x ] +
                                LATTICE[ y + n*PBC(x-1, n) ] +
                                LATTICE[ PBC(y-1, n) + n*x ]  ) * LATTICE[y + n*x];
        //field
        delE -= 2 * F * LATTICE[y + n*x];
        //if u_new less than initial u, flip the spin 
        if(delE < 0.) {
            LATTICE[y + n*x] *= -1;
        } 
        else if (exp(-beta * delE) >= rand_unif()) { 
            LATTICE[y + n*x] *= -1;
        }
    }
}

void mc_step_non_boltz(int* LATTICE, int n, float beta, float J, float F, int* cluster_map, int* tmpmap, int N, float w, float w_size ) { 
    int x, y = 0;
    int num_spins = n*n;
    float umbr_constr = 0.;
    int y1, y2, y3, y4 = 0;
    float H = 0.;
    bool near_nucleus, in_nucleus, at_surface, in_windo_when_incr, in_windo_when_decr;
    int c1, c2 = 0;
    float delE = 0;
    int old_pos=0;
    int lower_bound =  (int)(w - w_size/2);
    int upper_bound =  (int)(w + w_size/2);
    int N_id = get_largest_id(cluster_map, tmpmap, n);
    float M = 0.;
    float rnd = 0.;
    float E_UMBR = 0.;
    bool safe=false;
    srand(time(0));
    bool first = true;
    for(int i = 0; i < n*n; i++) { 
        x = (int) rand()%n;
        y = (int) rand()%n;
        y1 =  y + n*PBC(x+1, n);
        y2 =  PBC(y+1, n) + n*x;
        y3 =  y + n*PBC(x-1, n);
        y4 =  PBC(y-1, n) + n*x; 
        // main hamiltonian term
        delE =  2 * J * (LATTICE[y1]+LATTICE[y2]+LATTICE[y3]+LATTICE[y4]) * LATTICE[y + n*x];
        // field term
        delE += 2 * F * LATTICE[y + n*x];
        // constraints
        in_nucleus = (cluster_map[y + n*x] == N_id);
        near_nucleus = (cluster_map[y1] == N_id)||(cluster_map[y2] == N_id)||(cluster_map[y3] == N_id)||(cluster_map[y4] == N_id);
        //at_surface   = (cluster_map[y1] != N_id)||(cluster_map[y2] != N_id)||(cluster_map[y3] != N_id)||(cluster_map[y4] != N_id);
        in_windo_when_incr = (lower_bound < N + 1) && (upper_bound > N + 1);
        in_windo_when_decr = (lower_bound < N - 1) && (upper_bound > N - 1);
        // only remove if strongly connected to nucleus
        safe = check_edge(LATTICE, cluster_map, y+n*x, n, N_id);
        // metropolis mc step
        if ((delE < 0) || (exp(-beta * delE) > rand_unif())) {
            if(in_nucleus && in_windo_when_decr) { 
                LATTICE[y + n*x] *= -1; //spin flip operation
                cluster_map[y + n*x] = -1;
                N--;
            }
            else if(!in_nucleus && near_nucleus && in_windo_when_incr) { 
                LATTICE[y + n*x] *= -1; //spin flip operation
                cluster_map[y + n*x] = N_id;
                N++;
            }
        }
    }
}

