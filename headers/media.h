#ifndef MEDIA_H
#define MEDIA_H

#include <stdint.h>
#include <string.h>
#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <stdint.h>
#include <stdbool.h>
#include <time.h>
#include <string.h>
#include "../headers/util.h"

void save_snapshot(int *LATTICE, uint8_t* tmp,  int n, int i);
void save_video(int* simspace, int n, int nsteps);
void save_video_clr(int* simspace, int n, int nsteps, int* cluster_map, int* N_id);

#endif