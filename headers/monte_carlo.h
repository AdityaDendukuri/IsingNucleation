#ifndef MONTE_CARLO_H
#define MONTE_CARLO_H



void mc_step_boltz(int* LATTICE, int n,  float beta, float J, float F);
void mc_step_non_boltz(int* LATTICE, int n, float beta, float J, float F, int* cluster_map, int* tmpmap, int N, float w, float w_size );



#endif