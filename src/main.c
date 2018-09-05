#include <fftw3.h>
#include <hdf5.h>
#include <math.h>
#include <omp.h>
#include <stdio.h>
#include <stdlib.h>
#include <assert.h>
#include <time.h>

#include "allvars.h"
#include "utils.h"



fftw_plan plan_fft_psi;
fftw_plan plan_ifft_psi;
fftw_plan plan_fft_V;
fftw_plan plan_ifft_V;  


void kick(double dt);
void drift(double dt);
void updatePotential(void);
double get_kSq(int i, int j, int k);



/** Belerofon Main loop
 * 
 */
int main( int argc, char* argv[] ) {
  
  
    // set the number of threads from input argument (otherwise, use default of 4)
    Num_threads = (argc > 1) ? atoi(argv[1]) : 4; 
    
    omp_set_dynamic(0);     // Explicitly disable dynamic teams
    omp_set_num_threads(Num_threads);
  
    printf( "\n\nBELEROFON\n\n" );
    printf( "Starting simulation ...\n" );  
    printf("Using %d threads\n\n", Num_threads);


    // [0] init

    Lbox = 1.0 / BETA;
    a = get_a_from_snap(0);
    double t = get_t(a);
    make_sim_dir();
    printf("  Simulation start is at a=%0.2f, t=%0.2f\n", a, t);
    printf("  Simulation end is at a=%0.2f, t=%0.2f\n", get_a_from_snap(NOUT), get_t_from_snap(NOUT));

	int snap = 0;
	int saveNextTurn = 1;
	double dt = 0;
    readIC(); // updates the global variable a
    t = get_t(a);
    snap = set_snap(a);
    printf("IC is at a=%0.2f, t=%0.2f\n", a, t);
    printf("First snapshot to write is snap %d of %d\n\n", snap, NOUT);
    
    #pragma omp parallel for
	for(int i = 0; i < N; i++) {
		for(int j = 0; j < N; j++) {
			for(int k = 0; k < N; k++) {
				V[i][j][k] = 0; // initialize potential to 0
			}
		}
	}
	
	int fftw_init_threads(void);
	void fftw_plan_with_nthreads(int Num_threads);
	
	plan_fft_psi  = fftw_plan_dft_3d(N, N, N, psi, psihat, FFTW_FORWARD,  FFTW_ESTIMATE);
	plan_ifft_psi = fftw_plan_dft_3d(N, N, N, psihat, psi, FFTW_BACKWARD, FFTW_ESTIMATE);
	plan_fft_V    = fftw_plan_dft_r2c_3d(N, N, N, V, Vhat, FFTW_ESTIMATE);
	plan_ifft_V   = fftw_plan_dft_c2r_3d(N, N, N, Vhat, V, FFTW_ESTIMATE);
	
  
	// [1] Main loop
	
    clock_t tic, toc;
    tic = clock();
	
	double t_final = get_t_from_snap(NOUT);
	double t_nextsnap = get_t_from_snap(snap);
	int cc = 0;
	
	while (t < t_final) {
			
		// [kick-drift-(pot)-kick]
		kick(0.5*dt);		
		drift(dt);		
		updatePotential();		
		kick(0.5*dt);		

	
		t += dt;
		a = get_a(t);
		cc += 1;
		if ( cc % 100 == 0 ) {
			toc = clock();
			printf("    step: %d,   a=%0.3f,   took %0.2f s\n", cc, a, (double)(toc - tic) / CLOCKS_PER_SEC);
			tic = clock();
		}
		
	
		// save output
		if (saveNextTurn == 1) {
			saveSnap(snap);
			snap += 1;
			saveNextTurn = 0;
			t_nextsnap = get_t_from_snap(snap);
		}
	
		
		// set next timestep
		dt = get_dt(a);
		
		if (snap <= NOUT)
			if (t+dt >= t_nextsnap) {
				dt = t_nextsnap - t;
				assert(dt > 0);
				saveNextTurn = 1;
			}
			
	}
	// END main loop



	// clean up
	fftw_destroy_plan(plan_ifft_V);
	fftw_destroy_plan(plan_fft_V);
	fftw_destroy_plan(plan_ifft_psi);
	fftw_destroy_plan(plan_fft_psi);

	fftw_free(Vhat);
	fftw_free(V);
	fftw_free(psihat); 
	fftw_free(psi);
  
  
	printf( "Simulation completed!\n\n" );
		
	return 0;
}






void kick(double dt) {
	// [kick]    psi = exp(-1.i * dt/2 * V).*psi;
	#pragma omp parallel for
	for(int i = 0; i < N; i++) {
		for(int j = 0; j < N; j++) {
			for(int k = 0; k < N; k++) {
				double pRe = psi[i][j][k][0];
				double pIm = psi[i][j][k][1];
				double cosA = cos(-dt*V[i][j][k]);
				double sinA = sin(-dt*V[i][j][k]);
				psi[i][j][k][0] = cosA * pRe - sinA * pIm; // (c + I* s) * (pR + I* pI) == 
				psi[i][j][k][1] = sinA * pRe + cosA * pIm; 
			}
		}
	}	
}




void drift(double dt) {
	// [drift] psi = ifftn(exp(dt * (-1.i*kSq/2/a^2)).*fftn(psi));
	fftw_execute(plan_fft_psi);
    #pragma omp parallel for
	for(int i = 0; i < N; i++) {
		for(int j = 0; j < N; j++) {
			for(int k = 0; k < N; k++) {
				double kSq = get_kSq(i,j,k);
				double pRe = psihat[i][j][k][0];
				double pIm = psihat[i][j][k][1];
				double cosA = cos(-dt*0.5*kSq/a/a);
				double sinA = sin(-dt*0.5*kSq/a/a);
				psihat[i][j][k][0] = (cosA * pRe - sinA * pIm) / (N*N*N);
				psihat[i][j][k][1] = (sinA * pRe + cosA * pIm) / (N*N*N);
			}
		}
	}
	fftw_execute(plan_ifft_psi);
	
}




void updatePotential(void) {
	// [potential] V = ifftn(-fftn(G*beta^2/2*(abs(psi).^2 - 1)) ./ ( kSq  + (kSq==0))) / a - b * abs(psi).^2/a^3 ./ (1 + abs(psi).^2/a^3);
	
    #pragma omp parallel for
	for(int i = 0; i < N; i++) {
		for(int j = 0; j < N; j++) {
			for(int k = 0; k < N; k++) {
				V[i][j][k] = 0; // initialize potential to 0
				if (GSWITCH == 1) {
					double rho = psi[i][j][k][0] * psi[i][j][k][0] + psi[i][j][k][1] * psi[i][j][k][1];
					V[i][j][k] = 0.5*BETA*BETA*(rho - 1.0);
				}	
			}
		}
	}
	
	
	if (GSWITCH == 1) {
		fftw_execute(plan_fft_V);
		#pragma omp parallel for
		for(int i = 0; i < N; i++) {
			for(int j = 0; j < N; j++) {
				for(int k = 0; k < N/2+1; k++) {
					// after an r2c transform, the output is N x N x (N/2+1) 
					// i,j goes from 0,...,N-1 
					// k goes from 0,...,N/2   since Fhat(-k) = Fhat(k)^*
					double kSq = get_kSq(i,j,k);
					if (kSq == 0) kSq = 1.0;
					Vhat[i][j][k][0] = - Vhat[i][j][k][0] / kSq / a / (N*N*N);
					Vhat[i][j][k][1] = - Vhat[i][j][k][1] / kSq / a / (N*N*N);
				}
			}
		}
		fftw_execute(plan_ifft_V);
	}
	
	
	if (BSWITCH == 1) {
		#pragma omp parallel for
		for(int i = 0; i < N; i++) {
			for(int j = 0; j < N; j++) {
				for(int k = 0; k < N; k++) {
					double rho = psi[i][j][k][0] * psi[i][j][k][0] + psi[i][j][k][1] * psi[i][j][k][1];
					rho /= a*a*a;
					V[i][j][k] += - rho / (1.0 + rho);
				}
			}
		}
	}	
	
	
}




double get_kSq(int i, int j, int k) {
	// i,j,k goes from 0,...,N-1
	if (i >= N/2) i -= N;
	if (j >= N/2) j -= N;
	if (k >= N/2) k -= N;
	return ((double) i*i + j*j + k*k) * 4.0*M_PI*M_PI / Lbox / Lbox;
}


