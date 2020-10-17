#ifndef ISING_LATTICE_H
#define ISING_LATTICE_H

#include <stdbool.h>
#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <stdint.h>
#include <stdbool.h>
#include <time.h>
#include <limits.h>
#include <string.h>
#include "../headers/util.h"


float hamiltonian(int* LATTICE, int n, float J, float F);
void num_up(int* lattice, int n);
void random_init(int* LATTICE, int n);
int getid(int* LATTICE, int n);
float magnetization(int* LATTICE, int n);
float nucleation(int N, float del_mu, float gamma, float beta);
void  custom_init(int* LATTICE, int n);
int largest_cluster(int* adj, int n, int pos, int* cluster_map, int id);
int largest_cluster_size(int* LATTICE, int* adj, int* cluster_map, int n, int spin);
int get_largest_id(int* cluster_map, int* tmpmap, int n);
int surface_area(int* LATTICE, int* cluster_map, int n, int N_idx, int N);
bool check_centrality(int* LATTICE, int* cluster_map, int i, int n, int N_id);
bool check_safety(int* LATTICE, int* cluster_map, int i, int n, int N_id);
bool check_edge(int* LATTICE, int* cluster_map, int i, int n, int N_id);
int place_solvent(int* LATTICE, int n, int x, int y, int N, int n_flipped);
void generate_nucleus(int* LATTICE, int n, float r);

#endif
