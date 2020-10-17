#include "../headers/simul.h"
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


void run_sim_nucleation(int* LATTICE, int n,  float J, float F, float beta) {
    srand(time(0));
    int n_windows = 30;
    int num_spins = n*n;
    int num_init = 2000;
    int num_sweeps = 10000;
    int N_tot = 0;
    int N = 0;
    int N_id = 0;
    int N_window = 15;
    int window_size = 30;
    bool next_sample_chosen = false;
    int iter=0;
    int lower_bound = 0;
    int upper_bound = 0;
    float M = 0;
    int under = 0;
    int over = 0;
    int display_steps= 5;
    //define memory
    int* mag_space = malloc(num_sweeps * n_windows * sizeof(int));
    int* tmpmap = malloc(sizeof(int)*num_spins);
    int* adj = malloc(sizeof(int)*num_spins*4);
    int* cluster_map = malloc(sizeof(int) * num_spins);
    int* start_state = malloc(sizeof(int)*num_spins);
    int* original_state = malloc(sizeof(int)*num_spins);
    int* tmp_state = malloc(sizeof(int)*num_spins);
    int* simspace = malloc(sizeof(int)*num_sweeps*n_windows*num_spins/display_steps);
    //int* sim_space;
    //uint8_t* image = malloc(n*n);
    //assign initial values to memory 
    set(cluster_map, -1, num_spins);
    set(mag_space, -1, num_sweeps * n_windows);
    generate_nucleus(LATTICE, n, N_window + 500);
    memcpy(start_state, LATTICE, n*n*sizeof(int));
    memcpy(original_state, LATTICE, n*n*sizeof(int));
    int k = 0;
    int l = 0;
    int old_n = 0;
    //iterating each window
    for(int i = 0; i < n_windows; i++) { 
        lower_bound = (int)(N_window - window_size/2);
        upper_bound = (int)(N_window + window_size/2);
        printf( "WINDOW %d , N = %d: \n", i+1, N_window);
        printf("    Equilibrating : "); fflush(stdout);
        if(next_sample_chosen) { 
            memcpy(LATTICE, start_state, n*n*sizeof(int));
            next_sample_chosen = false;
        }
        else {
            printf("starting again ");
            memcpy(LATTICE, original_state, n*n*sizeof(int));
        }
        // equilibrate 
        for(int j = 0; j < num_init; j++) { 
            N = largest_cluster_size(LATTICE, adj, cluster_map, n, -1);
            if ( (N >= lower_bound)&&(N <= upper_bound)) break;
            mc_step_boltz(LATTICE, n, beta, J, 0);
        }
        //copy lattice state into two variables and calculate current N
        memcpy(start_state, LATTICE, n*n*sizeof(int));
        memcpy(original_state, LATTICE, n*n*sizeof(int));
        N = largest_cluster_size(LATTICE, adj, cluster_map, n, -1);
        M = magnetization(LATTICE, n);
        if ((N < lower_bound)||(N > upper_bound)) {
            break;
        }
        int j = 0;
        while(j < 2000) { 
            memcpy(tmp_state, LATTICE, n*n*sizeof(int));
            mc_step_boltz(LATTICE, n, beta, J, F);
            N = largest_cluster_size(LATTICE, adj, cluster_map, n, -1);
            if ( (N < lower_bound) || (N > upper_bound)) {
                memcpy(LATTICE, tmp_state, n*n*sizeof(int));
                N = old_n;
            } 
            old_n = N;
            if( j % 200 == 0) { 
                printf(" %d ", N); fflush(stdout);
            }
            j++;
        }
        N = largest_cluster_size(LATTICE, adj, cluster_map, n, -1);
         printf("\n");
        printf("NUC : %d \n", N);
        //save_snapshot(LATTICE, image, n, i);
        printf("    Simulation    : ");
        j = 0;
        over = 0;
        under = 0;

        while(j < num_sweeps) { 
            memcpy(tmp_state, LATTICE, n*n*sizeof(int));
            //mc_step_non_boltz(LATTICE,n, beta, J, F, cluster_map, tmpmap,N,N_window, window_size);
            mc_step_boltz(LATTICE, n, beta, J, F);
            //bool eq = equal(LATTICE, tmp_state, n);
            N = largest_cluster_size(LATTICE, adj, cluster_map, n, -1);
            if ( (N < lower_bound) || (N >  upper_bound)) {
                memcpy(LATTICE, tmp_state, n*n*sizeof(int));
                N = old_n;
            } 
            old_n = N;
            if(N > N_window) 
                over++;
            else 
                under++;
            mag_space[j + num_sweeps*i] = N;
            if((N >= N_window)&&(N <= upper_bound)) {
                memcpy(start_state, LATTICE, n*n*sizeof(int));
                next_sample_chosen = true;
            }
            if(k%display_steps==0) { 
                memcpy(&simspace[num_spins*l], LATTICE, num_spins*sizeof(int));
                l++;
            }
            j++;
            k++;
            if ( j % 1000 == 0) { 
                printf(" %d ", N); fflush(stdout);
            }
            
        }
        printf("\nBelow midpoint : %d percent \n",  under);
        printf("Above midpoint : %d percent \n", over);
        printf("\n \n \n");
        // move window 
        N_window += 15;
    }
    save_video(simspace, n, l);
    printf("VIDEO SAVED \n");
    FILE *fp = fopen("../output/mag.txt", "w+"); 
    for(int i = 0; i < num_sweeps * n_windows; i++) { 
        fprintf(fp, "%d \n", mag_space[i]);
    }
    printf("DATA SAVED \n");


    free(mag_space);
    free(cluster_map);
    free(start_state);
    free(tmpmap);
    free(tmp_state);
    //free(image);
    free(original_state);
    free(adj);
    free(simspace);
    fclose(fp);
    
}
void ones(int* LATTICE, int n) { 
    for(int i = 0; i < n*n; i++) { 
        LATTICE[i] =  1;
    }
}

void run_sim(int* LATTICE,  int n, float J, float F, float beta) { 
    int num_spins = n*n;
    FILE *fp_mag = fopen("../output/mag.txt", "w+");
    int num_sweeps = 1000;
    int init_steps = 0;
    int dstep = 100;
    int sstep = 1;
    int video_frequency = 1;
    float uinit = hamiltonian(LATTICE, n, J, F);
    float prob = 0.;
    float delE = 0.;
    float mag;
    int* simspace = malloc(sizeof(int)*num_sweeps*num_spins/video_frequency);
    int k = 0;
    random_init(LATTICE, n);
    // trial steps to get to equilibrium
    for(int i = 0; i < init_steps; i++) { 
        mc_step_boltz(LATTICE,n, beta, J, F);
    }
    // main simulation
    for(int i = 0; i < num_sweeps; i++) { 

        mc_step_boltz(LATTICE,n, beta, J, F);
        uinit = hamiltonian(LATTICE, n, J, F);
        mag = magnetization(LATTICE, n);
        //display attributes
        if(i%dstep == 0) {    
            num_up(LATTICE, n);
            printf("Sweep #%d -> U :%f | M:%f \n", i, uinit, mag);
	    }
        // save data 
        if(i%sstep==0) {  
            fprintf(fp_mag, "%f \n", mag);
        }
        if(i%video_frequency==0) { 
            memcpy(&simspace[num_spins*k], LATTICE, num_spins*sizeof(int));
            k++;
        }
    }
    //save simulation video
    save_video(simspace, n, k);
    // free all allocated memory and close all files 
    fclose(fp_mag);
}

void validate_free_energy(int* LATTICE, int n, int J, int F) { 
    int num_spins = n*n;
    int num_equib = 5000;
    int num_sweeps = 10000;
    int ssteps = 10; 
    int dsteps = 1000;
    int vid_frequency = 10;
    float beta = 0;
    float m = 0;
    float T_MIN = 1.2;
    float T_MAX = 2.4;
    int n_T = (int)ceil((T_MAX - T_MIN)/0.1);
    int* simspace = malloc(n*n*sizeof(int)*n_T*num_sweeps/vid_frequency);
    FILE *fp = fopen("../output/free_energy_test_result.txt", "w+");
    random_init(LATTICE, n);
    int k = 0;
    for(float T = T_MIN; T <= T_MAX; T+=0.1) { 
        beta = 1./T; 
        printf("T = %f : ", T);
        fflush(stdout);
        // equilibrating steps
        for(int i = 0; i < num_equib; i++) { 
            mc_step_boltz(LATTICE, n, beta, J, F);
        }
        // main sim loop
        for(int i = 0; i < num_sweeps; i++) { 
            mc_step_boltz(LATTICE, n, beta, J, F);
            if(i%ssteps==0) { 
                m = magnetization(LATTICE, n);
                fprintf(fp, "%f\n", fabs(m));
            }
            if(i%dsteps==0) {
                m = magnetization(LATTICE, n);
                printf(" %f ", m);
                fflush(stdout);
            }
            // add a frame to video every other step
            if(i % vid_frequency == 0) { 
                memcpy(&simspace[num_spins*k], LATTICE, num_spins*sizeof(int));
                k++;
            }
        }
        printf("\n");
    }
    save_video(simspace, n, k);
    free(simspace);
    fclose(fp);
}