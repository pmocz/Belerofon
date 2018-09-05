#include <fftw3.h>
#include <hdf5.h>
#include <math.h>
#include <omp.h>
#include <stdio.h>
#include <stdlib.h>
#include <assert.h>

#include "allvars.h"
#include "utils.h"



fftw_plan plan_fft_psi;
fftw_plan plan_ifft_psi;
fftw_plan plan_fft_V;
fftw_plan plan_ifft_V;  


void kick(double dt);
void drift(double dt);
void updatePotential(void);
double get_kSq(int idx);




/** 
 * BELEROFON
 * Philip Mocz
 * Princeton University (2018)
 * Gross-Pitaevskii-Poisson Solver
 */
int main( int argc, char* argv[] ) {
  
  
    // Set # threads from input argument (default to 4)
    Num_threads = (argc > 1) ? atoi(argv[1]) : 4; 
    
    omp_set_dynamic(0);     // Explicitly disable dynamic teams
    omp_set_num_threads(Num_threads);
  
    printf( "\n\nBELEROFON\n\n" );
    printf( "Gross-Pitaevskii-Poisson Solver\n\n\n" );
    printf( "Starting simulation ...\n" );  
    printf( "Using %d threads\n\n", Num_threads );


	int fftw_init_threads(void);
	void fftw_plan_with_nthreads(int Num_threads);

    // Allocate FFT arrays
	psi	   = (fftw_complex*) fftw_malloc(sizeof(fftw_complex) * N*N*N);
	psihat = (fftw_complex*) fftw_malloc(sizeof(fftw_complex) * N*N*N);
	V      = (double*)       fftw_malloc(sizeof(double)       * N*N*N);
	Vhat   = (fftw_complex*) fftw_malloc(sizeof(fftw_complex) * N*N*(N/2+1));
	
	plan_fft_psi  = fftw_plan_dft_3d(N, N, N, psi, psihat, FFTW_FORWARD,  FFTW_ESTIMATE);
	plan_ifft_psi = fftw_plan_dft_3d(N, N, N, psihat, psi, FFTW_BACKWARD, FFTW_ESTIMATE);
	plan_fft_V    = fftw_plan_dft_r2c_3d(N, N, N, V, Vhat, FFTW_ESTIMATE);
	plan_ifft_V   = fftw_plan_dft_c2r_3d(N, N, N, Vhat, V, FFTW_ESTIMATE);


    /* [0] Initialize */

    Lbox = 1.0 / BETA;
    
	make_sim_dir();
    ioSnap(READ_MODE,0); // read snapshot.  updates the global variable a
    
    double t = get_t(a);
    int snap = set_snap(a);
	int saveNextTurn = 1;
	double t_final = get_t_from_snap(NOUT);
	double t_nextsnap = get_t_from_snap(snap);
	int cc = 0;
	double dt = 0;
	double tic, toc;
	
    printf("  IC is at a=%0.2f, t=%0.2f\n", a, t);
    printf("  First snapshot to write is snap %d of %d\n\n", snap, NOUT);
    
	for(int i = 0; i < N*N*N; i++)
		V[i] = 0; // initialize potential to 0
	
  
	/* [1] Main Loop */ 
    tic = omp_get_wtime();
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
			toc = omp_get_wtime();
			printf("    step: %d,   a=%0.3f,   took %0.2f s\n", cc, a, (toc - tic));
			tic = omp_get_wtime();
		}
		
	
		// save output
		if (saveNextTurn == 1) {
			ioSnap(WRITE_MODE,snap);
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
	// END Main Loop
	toc = omp_get_wtime();


	// clean up
	fftw_destroy_plan(plan_ifft_V);
	fftw_destroy_plan(plan_fft_V);
	fftw_destroy_plan(plan_ifft_psi);
	fftw_destroy_plan(plan_fft_psi);
  
	printf( "Simulation completed!\n\n" );
		
	return 0;
}






void kick(double dt) {
	// [kick]    psi = exp(-1.i * dt/2 * V).*psi;
	int i;
	#pragma omp parallel for
	for(i = 0; i < N*N*N; i++) { 
		double pRe = psi[i][0];
		double pIm = psi[i][1];
		double cosA = cos(-dt*V[i]);
		double sinA = sin(-dt*V[i]);
		psi[i][0] = cosA * pRe - sinA * pIm; // (c + I* s) * (pR + I* pI) == 
		psi[i][1] = sinA * pRe + cosA * pIm; 
	}	
}




void drift(double dt) {
	// [drift] psi = ifftn(exp(dt * (-1.i*kSq/2/a^2)).*fftn(psi));
	fftw_execute(plan_fft_psi);
	int i;
	#pragma omp parallel for
	for(i = 0; i < N*N*N; i++) {
		double kSq = get_kSq(i);
		double pRe = psihat[i][0];
		double pIm = psihat[i][1];
		double cosA = cos(-dt*0.5*kSq/a/a);
		double sinA = sin(-dt*0.5*kSq/a/a);
		psihat[i][0] = (cosA * pRe - sinA * pIm) / (N*N*N);
		psihat[i][1] = (sinA * pRe + cosA * pIm) / (N*N*N);
	}
	fftw_execute(plan_ifft_psi);
	
}




void updatePotential(void) {
	// [potential] V = ifftn(-fftn(G*beta^2/2*(abs(psi).^2 - 1)) ./ ( kSq  + (kSq==0))) / a - b * abs(psi).^2/a^3 ./ (1 + abs(psi).^2/a^3);
	
	int i;
	#pragma omp parallel for
	for(i = 0; i < N*N*N; i++) {
		V[i] = 0; // initialize potential to 0
#if GSWITCH == 1
		double rho = psi[i][0] * psi[i][0] + psi[i][1] * psi[i][1];
		V[i] = 0.5*BETA*BETA*(rho - 1.0);
#endif	
	}
	
	
#if GSWITCH == 1
	fftw_execute(plan_fft_V);
	#pragma omp parallel for
	for(i = 0; i < N*N*(N/2+1); i++) {
		// after an r2c transform, the output is N x N x (N/2+1) 
		// i,j goes from 0,...,N-1 
		// k goes from 0,...,N/2   since Fhat(-k) = Fhat(k)^*
		double kSq = get_kSq(i);
		if (kSq == 0) kSq = 1.0;
		Vhat[i][0] = - Vhat[i][0] / kSq / a / (N*N*N);
		Vhat[i][1] = - Vhat[i][1] / kSq / a / (N*N*N);
	}
	fftw_execute(plan_ifft_V);
#endif
	
	
#if BSWITCH == 1
	#pragma omp parallel for
	for(i = 0; i < N*N*N; i++) {
		double rho = psi[i][0] * psi[i][0] + psi[i][1] * psi[i][1];
		rho /= a*a*a;
		V[i] += - rho / (1.0 + rho);
	}
#endif	
	
	
}




double get_kSq(int idx) {  // linear index idx = k + N * (j + N * i) = [i][j][k]
	int i = idx / (N*N);
	int j = (idx % (N*N)) / N;
	int k = idx % N;
	// i,j,k goes from 0,...,N-1
	if (i >= N/2) i -= N;
	if (j >= N/2) j -= N;
	if (k >= N/2) k -= N;
	return ((double) i*i + j*j + k*k) * 4.0*M_PI*M_PI / Lbox / Lbox;
}


