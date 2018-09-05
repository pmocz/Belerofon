#ifndef ALLVARS_H
#define ALLVARS_H

#include <fftw3.h>

/* Global Variables and Defines here */

// #####################################################################

/* PARAMETERS  */
#define N                64     // resolution              
#define BETA           0.06     // gravity    
#define BSWITCH           1     // self-interaction switch 
#define GSWITCH           1     // gravity switch
                                
#define NOUT             80     // number of output files  0..NOUT
#define LOGAFINAL       2.0     // simulate from a = 1..10^LOGAFINAL

#define ICFILE   "output/ic64box16.6667B0.06.h5"  // can restart from snapshot too!


// #####################################################################


#define KMAX     0.70710678118  // sqrt(2)/2
#define SQRT3P2  0.86602540378  // sqrt(3)/2
#define M_PI     3.14159265358979323846
#define READ_MODE     0
#define WRITE_MODE    1

int Num_threads;

double Lbox;

char SimDir[128];

double a;

fftw_complex *psi;
fftw_complex *psihat;
double *V;
fftw_complex *Vhat;

#endif
