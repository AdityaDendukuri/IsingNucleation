#include "../headers/util.h"



int PBC(int a, int b) {
  return (((a % b) + b) % b);
}

float sigmoid(float x) { 
    return 1. / (1. + exp(-x));
}

void set(int* x, int val, int n) { 
    for(int i = 0; i < n; i++) { 
        x[i] = val;
    }
}


double rand_unif() {
   // return a uniformly distributed random value
   return ((double)(rand()))/( (double)(RAND_MAX) + 1. );
}

double normalRandom() {
   // return a normally distributed random value
   double v1=rand_unif();
   double v2=rand_unif();
   return sigmoid(cos(2*3.14*v2)*sqrt(-2.*log(v1)));
}



