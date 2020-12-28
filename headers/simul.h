#ifndef SIMUL_H
#define SIMUL_H



void run_sim_nucleation(int* LATTICE, int n,  float J, float F, float beta);
void ones(int* LATTICE, int n);
void run_sim(int* LATTICE,  int n, float J, float F, float beta);
void validate_free_energy(int* LATTICE, int n, int J, int F);
void make_nucleus(int* LATTICE,  int n, float J, float F, float beta);
#endif