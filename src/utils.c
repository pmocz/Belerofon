#include <sys/types.h>
#include <sys/stat.h>
#include <unistd.h>
#include <stdio.h>
#include <string.h>
#include <stdlib.h>
#include <math.h>
#include <fftw3.h>
#include <hdf5.h>
#include <assert.h>
#include <omp.h>

#include "utils.h"
#include "allvars.h"


void make_sim_dir(void) {

	struct stat st = {0};
	char buffer[8];

    strncpy(SimDir, "output/N", sizeof(SimDir));
    sprintf(buffer,"%d",N);
    strncat(SimDir, buffer, (sizeof(SimDir) - strlen(SimDir)) );
    strncat(SimDir, "box", (sizeof(SimDir) - strlen(SimDir)) );
    sprintf(buffer,"%0.4f",Lbox);
    strncat(SimDir, buffer, (sizeof(SimDir) - strlen(SimDir)) );
    strncat(SimDir, "B", (sizeof(SimDir) - strlen(SimDir)) );
    sprintf(buffer,"%0.2f",BETA);
    strncat(SimDir, buffer, (sizeof(SimDir) - strlen(SimDir)) );
    strncat(SimDir, "b", (sizeof(SimDir) - strlen(SimDir)) );
    sprintf(buffer,"%d",BSWITCH);
    strncat(SimDir, buffer, (sizeof(SimDir) - strlen(SimDir)) );
    strncat(SimDir, "G", (sizeof(SimDir) - strlen(SimDir)) );
    sprintf(buffer,"%d",GSWITCH);
    strncat(SimDir, buffer, (sizeof(SimDir) - strlen(SimDir)) );
    strncat(SimDir, "/", (sizeof(SimDir) - strlen(SimDir)) );

    printf("output directory: %s\n", SimDir);

	if (stat("output/", &st) == -1)
		mkdir("output/", 0700);

	if (stat(SimDir, &st) == -1)
		mkdir(SimDir, 0700);

}


double get_a_from_snap(int snap) {
	return pow(10., (double) snap * LOGAFINAL / NOUT );
}

double get_t_from_snap(int snap) {
	return get_t( get_a_from_snap(snap) );
}

double get_t(double a) {
	return pow(a, 1.5) / (SQRT3P2*BETA);
}

double get_a(double t) {
	return pow( SQRT3P2*BETA*t, 0.66666666667 );
}

int set_snap(double a) {
	return round( log10(a) * ((double) NOUT / LOGAFINAL) );
}

double get_dt(double a) {
	double dx = Lbox / N;
	double dt = dx * dx / 6.0 * a * a;
	double Vmax = 1.0e-16;
	
	// single thread version
	//for(int i = 0; i < N*N*N; i++) 
	//	if (fabs(V[i]) > Vmax) Vmax = fabs(V[i]);

	double Vmax_thread[Num_threads];
	
	#pragma omp parallel
	{
		int id, n, start, stop;
		double my_max;
		
		n = N*N*N/Num_threads;         // step = N/number of threads.
		id = omp_get_thread_num(); // id is one of 0,1,...,(Num_threads-1)
		
		start = id * n;
		stop  = ( id != (Num_threads-1) )  ?  start + n : N;
		
		my_max = V[start];
		
		for (int i = start+1; i < stop; i++ )
			if ( fabs(V[i]) > my_max )
				my_max = fabs(V[i]) ;
		
		Vmax_thread[id] = my_max;	// Store result in min[id]   
	}
      
	Vmax = Vmax_thread[0];

	for (int i = 1; i < Num_threads; i++)
		if ( Vmax_thread[i] > Vmax )
			Vmax = Vmax_thread[i];
      
	
	dt = fmin(dt, 1.0 / Vmax);
	return dt;
}



void ioSnap(int mode, int snap) {
	// set filename
	struct stat st = {0};
	char buffer[8];
	char fname[128];
	
	// hdf5 handles
    hid_t    file_id, dataset_psiRe, dataset_psiIm, dataset_a;
    hid_t    dataspace_a, dataspace_psiRe, dataspace_psiIm;
    hid_t    memspace; 
    herr_t   status;
    
    // memspace variables
    hsize_t  mdim[2]         = {N*N*N,2}; // memory space dimensions
    hsize_t  count_out[2]    = {N*N*N,1}; // size of the hyperslab in memory
    hsize_t  stride_out[2]   = {1,2};     // hyperslab stride in memory
    hsize_t  offset_outRe[2] = {0,0};     // hyperslab offset in memory
    hsize_t  offset_outIm[2] = {0,1};     // hyperslab offset in memory

	// dataset variables
	hsize_t adim[1] = {1};
	hsize_t dims[3] = {N, N, N};
	
	
	// load/create hdf5 file
	if(mode == READ_MODE) {
		
		printf( "Reading IC ...\n" );
		
		file_id = H5Fopen(ICFILE, H5F_ACC_RDONLY, H5P_DEFAULT);
		
	} else {
		
		printf("Saving snap %d...\n", snap);
		strncpy(fname, SimDir, sizeof(fname));
		sprintf(buffer,"snap_");
		strncat(fname, buffer, (sizeof(fname) - strlen(SimDir)) );
		sprintf(buffer,"%03d",snap);
		strncat(fname, buffer, (sizeof(fname) - strlen(SimDir)) );
		sprintf(buffer,".h5");
		strncat(fname, buffer, (sizeof(fname) - strlen(SimDir)) );
		printf("  filename= %s\n", fname);
		
		if (stat(fname, &st) == 0)
			remove(fname);
		file_id = H5Fcreate(fname, H5F_ACC_TRUNC, H5P_DEFAULT, H5P_DEFAULT);
		
	}
	


    // [0] read/save a
    if(mode == READ_MODE) {
		
		dataset_a = H5Dopen(file_id, "/a", H5P_DEFAULT);
		status = H5Dread(dataset_a, H5T_NATIVE_DOUBLE, H5S_ALL, H5S_ALL, H5P_DEFAULT, &a);
		H5Dclose(dataset_a);
		
	} else {
		
		dataspace_a = H5Screate_simple(1, adim, NULL);
		dataset_a = H5Dcreate(file_id, "/a", H5T_NATIVE_DOUBLE, dataspace_a,
				H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT);
		status = H5Dwrite(dataset_a, H5T_NATIVE_DOUBLE, H5S_ALL, H5S_ALL,
				H5P_DEFAULT, &a);
		H5Dclose (dataset_a);
		H5Sclose (dataspace_a);
		
	}
	
	assert(status != -1);
    
    // [1] read/save psiRe
	memspace = H5Screate_simple(2, mdim, NULL);    // (Define memory hyperslab)
	status = H5Sselect_hyperslab(memspace, H5S_SELECT_SET, offset_outRe, stride_out, 
					count_out, NULL);
    
    if(mode == READ_MODE) {
		
		dataset_psiRe = H5Dopen(file_id, "/psiRe", H5P_DEFAULT);

		status = H5Dread(dataset_psiRe, H5T_NATIVE_DOUBLE, memspace, H5S_ALL,
				H5P_DEFAULT, psi);
	
		H5Dclose(dataset_psiRe);
		
	} else {

		dataspace_psiRe = H5Screate_simple(3, dims, NULL);
		dataset_psiRe = H5Dcreate(file_id, "/psiRe", H5T_NATIVE_DOUBLE, dataspace_psiRe,
				H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT);
				
		status = H5Dwrite(dataset_psiRe, H5T_NATIVE_DOUBLE, memspace, H5S_ALL,
				H5P_DEFAULT, psi);
		
		H5Dclose (dataset_psiRe);
		H5Sclose (dataspace_psiRe);
		
	}
	H5Sclose(memspace);
    
    // [2] read/save psiIm
	memspace = H5Screate_simple(2, mdim, NULL);    // (Define memory hyperslab)
	status = H5Sselect_hyperslab(memspace, H5S_SELECT_SET, offset_outIm, stride_out, 
					count_out, NULL);
    
    if(mode == READ_MODE) {
		
		dataset_psiIm = H5Dopen(file_id, "/psiIm", H5P_DEFAULT);

		status = H5Dread(dataset_psiIm, H5T_NATIVE_DOUBLE, memspace, H5S_ALL,
				H5P_DEFAULT, psi);
	
		H5Dclose(dataset_psiIm);
		
	} else {

		dataspace_psiIm = H5Screate_simple(3, dims, NULL);
		dataset_psiIm = H5Dcreate(file_id, "/psiIm", H5T_NATIVE_DOUBLE, dataspace_psiIm,
				H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT);
				
		status = H5Dwrite(dataset_psiIm, H5T_NATIVE_DOUBLE, memspace, H5S_ALL,
				H5P_DEFAULT, psi);
		
		H5Dclose (dataset_psiIm);
		H5Sclose (dataspace_psiIm);
		
	}
	H5Sclose(memspace);

    // [3] close file
    H5Fclose (file_id);

	printf( "i/o complete.\n" );	
}


