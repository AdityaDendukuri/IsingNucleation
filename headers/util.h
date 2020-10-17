#ifndef UTIL_H
#define UTIL_H

#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <stdint.h>
#include <stdbool.h>
#include <time.h>
#include <limits.h>
#include <string.h>
#include <string.h>

typedef unsigned char uchar;
int PBC(int a, int b);
float sigmoid(float x);
void set(int* x, int val, int n);
double rand_unif();
double normalRandom();

#endif