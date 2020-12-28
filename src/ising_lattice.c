#include "../headers/ising_lattice.h"


// Ising Spin Model hamiltonian
float hamiltonian(int* LATTICE, int n, float J, float F) {
    float H = 0.;
    int xidx = 0;
    int yidx = 0;
    int num_spins = n*n;
    for(int i = 0; i < num_spins; i++) {
        xidx = (int) i/n;
        yidx = (int) i%n;
        H += F * (LATTICE[i]);
        // access neighbours
        H += J * (  LATTICE[ yidx + n*PBC(xidx+1, n) ] +
                    LATTICE[ PBC(yidx+1, n) + n*xidx ] +
                    LATTICE[ yidx + n*PBC(xidx-1, n) ] +
                    LATTICE[ PBC(yidx-1, n) + n*xidx ]   ) * LATTICE[i];   
    }
    return -H/(4.*num_spins); 
}


void num_up(int* lattice, int n) { 
    int up = 0;
    int down = 0;
    for(int i = 0; i < n*n; i++) { 
        if (lattice[i] == 1) {up++;}
        if (lattice[i] == -1) {down++;}
    }
    printf("up : %d, down : %d \n ", up, down);
}



void random_init(int* LATTICE, int n) { 
    srand(time(0));
    float random=0.;
    for(int i = 0; i < n*n; i++) { 
        random = rand()%2;
        if(random==0){ LATTICE[i] = -1;}
        else{LATTICE[i] = +1;}
    }
}


// get structure id for integer
int getid(int* LATTICE, int n) { 
    int x = 0; 
    for(int i = 0; i < n; i++) { 
        if(LATTICE[i] != -1)
            x += pow(2, i);
    }
    return x;
}

float magnetization(int* LATTICE, int n) { 
    float M = 0.;
    int n_spins = n*n;
    for(int i = 0; i <n_spins; i++) {
        M += (float)LATTICE[i];
    }
    return M/n_spins;
}

float nucleation(int N, float del_mu, float gamma, float beta) { 
    return (- N * fabs(del_mu) + pow(N, 1/3)*gamma)*(1./beta);
}

void  custom_init(int* LATTICE, int n) { 
    memset(LATTICE, 1, n*n*sizeof(int));
    for(int i = 30; i < 90; i++) { 
        LATTICE[i + 50*n] = -1;
    }
    for(int i = 50; i < 90; i++) { 
        LATTICE[i + 30*n] = -1;
    }
    for(int i = 30; i < 100; i++) { 
        LATTICE[i + 90*n] = -1;
    }
}


int largest_cluster(int* adj, int n, int pos, int* cluster_map, int id) { 
    int n_clstr = 1;
    cluster_map[pos] = id;

    int x = pos/n;
    int y = pos%n;
    int y1 = y + n*PBC(x+1, n);
    int y2 = PBC(y+1, n) + n*x;
    int y3 = y + n*PBC(x-1, n);
    int y4 = PBC(y-1, n) + n*x;

    if( (cluster_map[y1] == -1) && (adj[0 + 4*pos] == 0) ) 
        n_clstr += largest_cluster(adj, n, y1, cluster_map, id);
    if( (cluster_map[y2] == -1) && (adj[1 +  4*pos] == 0) )
        n_clstr += largest_cluster(adj, n, y2, cluster_map, id);
    if( (cluster_map[y3] == -1) && (adj[2 +  4*pos] == 0) )
        n_clstr += largest_cluster(adj, n, y3, cluster_map, id);
    if( (cluster_map[y4] == -1) && (adj[3 +  4*pos] == 0) )
        n_clstr += largest_cluster(adj, n, y4, cluster_map, id);

    return n_clstr;
}


int largest_cluster_size(int* LATTICE, int* adj, int* cluster_map, int n, int spin) {
    int num_spins = n*n;
    int y1, y2, y3, y4 = 0;
    int x, y = 0;
    int cluster_size_tmp = 0;
    int cluster_size=0;
    int n_visited = 0;
    int id = 0;
    // adjacency matrix 
    memset(cluster_map, -1, sizeof(int)*num_spins);
    memset(adj, 1, sizeof(int)*num_spins*4);
    for(int i = 0; i < num_spins; i++) { 
        x = i/n;
        y = i%n;
        y1 = y + n*PBC(x+1, n);
        y2 = PBC(y+1, n) + n*x;
        y3 = y + n*PBC(x-1, n);
        y4 = PBC(y-1, n) + n*x;
        adj[0 + i*4] = fabs(LATTICE[i] - LATTICE[y1])/2;
        adj[1 + i*4] = fabs(LATTICE[i] - LATTICE[y2])/2;
        adj[2 + i*4] = fabs(LATTICE[i] - LATTICE[y3])/2;
        adj[3 + i*4] = fabs(LATTICE[i] - LATTICE[y4])/2;
    }
    //display_lattice(LATTICE, n);
    for(int i = 0; i < num_spins; i++) { 
        cluster_size_tmp = 0;
        if ((cluster_map[i] == -1)&&(LATTICE[i] == spin)) { 
            cluster_size_tmp = largest_cluster(adj, n, i, cluster_map, id);
            id++;
        }
        if(cluster_size < cluster_size_tmp) { 
            cluster_size = cluster_size_tmp;
        } 
    }
    return cluster_size;
}   
             /*             */
int get_largest_id(int* cluster_map, int* tmpmap, int n) { 
    memset(tmpmap, 0, sizeof(int)*n*n);
    for(int i = 0; i < n*n; i++) { 
        if(cluster_map[i] != -1)
            tmpmap[cluster_map[i]] += 1;
    }
    int largest_id = 0;
    for(int i = 0; i < n*n; i++) { 
        if(tmpmap[i] > tmpmap[largest_id]) { 
            largest_id = i;
        }
    }
    return largest_id;
}


int surface_area(int* LATTICE, int* cluster_map, int n, int N_idx, int N) { 
    int area = 0;
    int x, y, y1, y2, y3 , y4 = 0;
    bool at_surface;
    for(int i = 0; i < n*n; i++) { 
        area = 0;
        x  = i%n;
        y  = i%n;
        y1 =  y + n*PBC(x+1, n);
        y2 =  PBC(y+1, n) + n*x;
        y3 =  y + n*PBC(x-1, n);
        y4 =  PBC(y-1, n) + n*x; 
        at_surface = !((LATTICE[y1] == LATTICE[y2]) && (LATTICE[y2] == LATTICE[y3]) && (LATTICE[y3] == LATTICE[y4]));
        if ( (cluster_map[i] == N_idx) && (!at_surface) ) { 
            area += 1;
        }
    }
    return N - area;
}





bool check_centrality(int* LATTICE, int* cluster_map, int i, int n, int N_id) { 
    int x = (int) i/n;
    int y = (int) i%n;
    int y1 =  y + n*PBC(x+1, n);
    int y2 =  PBC(y+1, n) + n*x;
    int y3 =  y + n*PBC(x-1, n);
    int y4 =  PBC(y-1, n) + n*x;
    bool in_cluster = (cluster_map[i] == N_id);
    bool connected = (cluster_map[y1]==N_id)||(cluster_map[y2]==N_id)||(cluster_map[y3]==N_id)||(cluster_map[y4]==N_id);
    return (in_cluster&&connected);
}

bool check_safety(int* LATTICE, int* cluster_map, int i, int n, int N_id) { 
    cluster_map[i] = -1;
    int x = (int) i/n;
    int y = (int) i%n;
    int y1 =  y + n*PBC(x+1, n);
    int y2 =  PBC(y+1, n) + n*x;
    int y3 =  y + n*PBC(x-1, n);
    int y4 =  PBC(y-1, n) + n*x;
    bool check1 = check_centrality(LATTICE, cluster_map, y1, n, N_id);
    bool check2 = check_centrality(LATTICE, cluster_map, y2, n, N_id);
    bool check3 = check_centrality(LATTICE, cluster_map, y3, n, N_id);
    bool check4 = check_centrality(LATTICE, cluster_map, y4, n, N_id);
    cluster_map[i] = N_id;
    return check1 && check2 && check3 && check4;
}

bool check_edge(int* LATTICE, int* cluster_map, int i, int n, int N_id) { 
    int x = (int) i/n;
    int y = (int) i%n;
    int y1 =  PBC(y+1, n) + n*x;
    int y2 =  y + n*PBC(x+1, n);
    int y3 =  PBC(y-1, n) + n*x;
    int y4 =  y + n*PBC(x-1, n);
    bool check1 = ((cluster_map[y1] == N_id)&&(cluster_map[y2] == N_id));
    bool check2 = ((cluster_map[y2] == N_id)&&(cluster_map[y3] == N_id));
    bool check3 = ((cluster_map[y3] == N_id)&&(cluster_map[y4] == N_id));
    bool check4 = ((cluster_map[y4] == N_id)&&(cluster_map[y1] == N_id));
    return check1 || check2 || check3 || check4;

}




int place_solvent(int* LATTICE, int n, int x, int y, int N, int n_flipped) {
    if(n_flipped == N) {
        return n_flipped;
    }
    LATTICE[y + n*x] = -1;
    n_flipped++;
    int rnd = rand()%4;
    if(rnd==0) { 
        place_solvent( LATTICE, n,  PBC(x+1, n), y, N, n_flipped);
    }
    else if(rnd==1) { 
        place_solvent( LATTICE, n,  x, PBC(y+1, n), N, n_flipped);
    }
    else if(rnd==2) { 
        place_solvent( LATTICE, n,  PBC(x-1, n), y, N, n_flipped);
    }
    else if(rnd==3) { 
        place_solvent( LATTICE, n,  x, PBC(y-1, n), N, n_flipped);
    }
    return n_flipped;
}

void generate_nucleus(int* LATTICE, int n, float r) { 
    int center = n/2;
    int x = 0;
    int y = 0;
    float dist = 0;
    int start = (int)n/2-(int)sqrt(r)/2;
    int end = (int)n/2+(int)sqrt(r)/2;
    //memset(LATTICE, 1, sizeof(int)*n*n);
    for(int i = 0; i < n*n; i++) { 
        LATTICE[i] = 1;
    }
    for(int i = start; i < end; i++) { 
        for(int j = start; j < end; j++) { 
            LATTICE[j + n*i] = -1;
        }
    }
}