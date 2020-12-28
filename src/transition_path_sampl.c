#include "../headers/transition_path_sampl.h"
#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <stdint.h>
#include <stdbool.h>
#include <time.h>
#include <limits.h>
#include <string.h>
#include "../headers/ising_lattice.h"
#include "../headers/monte_carlo.h"
#include "../headers/media.h"
#include "../headers/util.h"


void trans_path_sampl(int* LATTICE, int n, int J, int F, float beta) { 
    srand(time(0));
    int n_windows = 30;
    int num_spins = n*n;
    int num_init = 25000;
    int lower_bound = 70;
    int upper_bound = 170;
    int N_tot = 0;
    int N = 0, M=0;
    int N_id = 0;
    int n_traj = 80;
    int traj_sweeps = 150;
    int disp_sweeps = floor(traj_sweeps / 10);
    int random = 0;
    //counters for reaction and product regions
    int num_selected = 0;
    int reac = 0;
    int prod = 0;
    //memory bufferes
    int* initial_states = malloc(num_init*num_spins*sizeof(int));
    int* state_picked = malloc(num_init*sizeof(int));
    int* tran_state_ensmb = malloc(n_traj*traj_sweeps*num_spins*sizeof(int)); 
    int* tran_state_ensmb_mp = malloc(n_traj*traj_sweeps*num_spins*sizeof(int)); 
    int* tran_state_C = malloc(n_traj*traj_sweeps*sizeof(int)); 
    int* tran_state_M = malloc(n_traj*traj_sweeps*sizeof(int));
    int* tran_state_N = malloc(n_traj*traj_sweeps*sizeof(int)); 
    int* adj = malloc(num_spins*4*sizeof(int));
    int* cluster_map = malloc(num_spins*sizeof(int));
    int* tmpmap = malloc(num_spins*sizeof(int));
    int* tmpstate = malloc(num_spins*sizeof(int));
    //initial steps to generate a small nucleus
    int k = 0;
    int old_n = 0;
    int Nid=0;
    set(state_picked, num_init, 0);
    generate_nucleus(LATTICE, n, 500);
    for(int j = 0; j < num_init; j++) { 
        N = largest_cluster_size(LATTICE, adj, cluster_map, n, -1);
        if ( (N >= lower_bound)&&(N <= upper_bound)) break;
        mc_step_boltz(LATTICE, n, beta, J, -0.08);
    }
    printf("INITIAL SWEEPS (step/N) : ");
    fflush(stdout);
    //initial steps
    for(int j = 0; j < num_init; j++) { 
        memcpy(tmpstate, LATTICE, n*n*sizeof(int));
        mc_step_boltz(LATTICE, n, beta, J, F);
        N = largest_cluster_size(LATTICE, adj, cluster_map, n, -1);
        if ( (N < lower_bound) || (N > upper_bound)) {
            memcpy(LATTICE, tmpstate, n*n*sizeof(int));
            N = old_n;
        } 
        old_n = N;
        if( j % (int)(num_init/10) == 0) { 
            printf(" %d ", N); fflush(stdout);
        }
        memcpy(&initial_states[j*num_spins], LATTICE, num_spins*sizeof(int));
    }
    //transition path ssampling loop
    printf("    BEGINNING TSP  with init N : %d \n", N);
    while(num_selected < n_traj) {
        //printf("\n ----------------------------------- \n");
        // pick a random time in the initial
        random = rand()%num_init;
        if(num_selected%10==0) { 
            //printf("\n REACTANT  : %d, PRODUCT : %d", reac, prod);
            printf("%d\n", num_selected);
        }
        state_picked[random] = 1;
        //save initial state 
        memcpy(LATTICE, &initial_states[random*num_spins], num_spins*sizeof(int));
        memcpy(tmpstate, LATTICE, num_spins*sizeof(int));
        //printf("STATE PICKED N %d : ", largest_cluster_size(LATTICE, adj, cluster_map, n, -1));
        //reset reactant/product count to zero
        reac = 0;
        prod = 0;
        //printf("Traj (sweep/N) : ");
        fflush(stdout);
        for(int i = 0; i < traj_sweeps; i++) { 
            //two forward mc steps is one shooting step
            mc_step_boltz(LATTICE, n, beta, J, F);
            //calculate largest cluster size  
            N = largest_cluster_size(LATTICE, adj, cluster_map, n, -1);
            M = magnetization(LATTICE, n);
            if(N < lower_bound) { 
                reac++;
                //memcpy(LATTICE, tmpstate, num_spins*sizeof(int));
            }
            else if(N > upper_bound) { 
                prod++;
                //memcpy(LATTICE, tmpstate, num_spins*sizeof(int));
            }
            //display
            //if( (i/2) % disp_sweeps == 0) { 
            //    printf(" %d/%d ", i, N);
            //    fflush(stdout);
            //}
            tran_state_N[k] = N;
            tran_state_M[k] = M;
            tran_state_C[k] = get_largest_id(cluster_map, tmpmap, n);
            memcpy(&tran_state_ensmb_mp[k*num_spins], cluster_map, num_spins*sizeof(int));
            memcpy(&tran_state_ensmb[k*num_spins], LATTICE, num_spins*sizeof(int));
            k++;
        }

        //if(num_selected%100 == 0)
        if(  reac > 0 && prod > 0 ) { 
            num_selected++;
            if(num_selected%10==0) { 
                printf("\n REACTANT  : %d, PRODUCT : %d \n", reac, prod);
                printf("%d\n\n", num_selected);
            }
            //printf(" ** TRAJECTORY ACCEPTED ** \n");
        }
        else { 
            //memcpy(LATTICE, tmpstate, num_spins*sizeof(int));
            k-=traj_sweeps;
            //printf(" ** TRAJECTORY REJECTED ** \n");
        }

        

        //printf("\n \n");
    }
    //save data
    FILE* fM = fopen("../output/tps_M.txt", "w+");
    FILE* fN = fopen("../output/tps_N.txt", "w+");
    for(int i = 0; i < n_traj*traj_sweeps; i++) { 
        fprintf(fM, "%d\n", tran_state_M[i]);
        fprintf(fN, "%d\n", tran_state_N[i]);
    }
    fclose(fN);
    fclose(fM);    
    printf("\n DATA SAVED :)\n");
    //save video
    save_video_clr(tran_state_ensmb, n, k, tran_state_ensmb_mp, tran_state_C);
    printf("\n Video SAVED :)\n");
    //free all mem
    free(initial_states);
    free(tran_state_ensmb);
    free(tran_state_ensmb_mp);
    free(tran_state_M);
    free(tran_state_C);
    free(tran_state_N);
    free(tmpmap);
    free(adj);
    free(cluster_map);
}