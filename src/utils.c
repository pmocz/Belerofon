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
	return pow(10.,(double) snap * LOGAFINAL / NOUT );
}

double get_t_from_snap(int snap) {
	return get_t(get_a_from_snap(snap));
}

double get_t(double a) {
	return pow(a, 1.5) / (SQRT3P2*BETA);
}

double get_a(double t) {
	return pow( SQRT3P2*BETA*t, 0.66666666667 );
}

double get_dt(double a) {
	double dx = Lbox / N;
	double dt = dx * dx / 6.0 * a * a;
	double Vmax = 1.0e-16;
	
	// single thread version
	//for(int i = 0; i < N; i++) {
	//	for(int j = 0; j < N; j++) {
	//		for(int k = 0; k < N; k++) {
	//          if (fabs(V[i][j][k]) > Vmax) Vmax = fabs(V[i][j][k]);
	//		}
	//	}
	//}
	
	double Vmax_thread[Num_threads];
	
	#pragma omp parallel
	{
		int id, n, start, stop;
		double my_max;
		
		n = N/Num_threads;         // step = N/number of threads.
		id = omp_get_thread_num(); // id is one of 0,1,...,(Num_threads-1)
		
		start = id * n;
		stop  = ( id != (Num_threads-1) )  ?  start + n : N;
		
		my_max = V[start][0][0];
		
		for (int i = start+1; i < stop; i++ )
			for(int j = 0; j < N; j++)
				for(int k = 0; k < N; k++) 
					if ( fabs(V[i][j][k]) > my_max )
						my_max = fabs(V[i][j][k]) ;
		
		Vmax_thread[id] = my_max;	// Store result in min[id]   
	}
      
	Vmax = Vmax_thread[0];

	for (int i = 1; i < Num_threads; i++)
		if ( Vmax_thread[i] > Vmax )
			Vmax = Vmax_thread[i];
      
	
	dt = fmin(dt, 1.0 / Vmax);
	return dt;
}


void readIC(void) {
	printf( "Reading IC ...\n" );
	
	// hdf5 handles
    hid_t       file_id, dataset_psiRe, dataset_psiIm, dataset_a; // handles
    hid_t       memspace; 
    herr_t      status;
    
    // memspace variables
    hsize_t  dimsm[4] = {N,N,N,2};        // memory space dimensions
    hsize_t  count_out[4]    = {N,N,N,1}; // size of the hyperslab in memory
    hsize_t  stride_out[4]   = {1,1,1,2}; // hyperslab stride in memory
    hsize_t  offset_outRe[4] = {0,0,0,0}; // hyperslab offset in memory
    hsize_t  offset_outIm[4] = {0,0,0,1}; // hyperslab offset in memory

 
    // Open the file and the datasets
    file_id = H5Fopen(ICFILE, H5F_ACC_RDONLY, H5P_DEFAULT);
    
    // [0] read a
    dataset_a = H5Dopen(file_id, "/a", H5P_DEFAULT);
    status = H5Dread(dataset_a, H5T_NATIVE_DOUBLE, H5S_ALL, H5S_ALL, H5P_DEFAULT, &a);
    assert(status != -1);
    H5Dclose(dataset_a);
    
    // [1] read psiRe
    dataset_psiRe = H5Dopen(file_id, "/psiRe", H5P_DEFAULT);
 
    // (Define memory hyperslab)
    memspace = H5Screate_simple(4,dimsm,NULL);   
    status = H5Sselect_hyperslab(memspace, H5S_SELECT_SET, offset_outRe, stride_out, 
				 count_out, NULL);

    status = H5Dread(dataset_psiRe, H5T_NATIVE_DOUBLE, memspace, H5S_ALL,
		     H5P_DEFAULT, psi[0]);

    H5Sclose(memspace);
    H5Dclose(dataset_psiRe);
 
    
    // [2] read psiIm
    dataset_psiIm = H5Dopen(file_id, "/psiIm", H5P_DEFAULT);
 
    memspace = H5Screate_simple(4,dimsm,NULL);   
    status = H5Sselect_hyperslab(memspace, H5S_SELECT_SET, offset_outIm, stride_out, 
				 count_out, NULL);

    status = H5Dread(dataset_psiIm, H5T_NATIVE_DOUBLE, memspace, H5S_ALL,
		     H5P_DEFAULT, psi[0]);

    H5Sclose(memspace);
    H5Dclose(dataset_psiIm);
    
    
    // close file
    H5Fclose(file_id);
	
	printf( "IC read.\n\n" );
}

void saveSnap(int snap) {
	printf("Saving snap %d...\n", snap);
	
	// set filename
	struct stat st = {0};
	char buffer[8];
	char fname[128];
	
	strncpy(fname, SimDir, sizeof(fname));
	sprintf(buffer,"snap_");
    strncat(fname, buffer, (sizeof(fname) - strlen(SimDir)) );
    sprintf(buffer,"%03d",snap);
    strncat(fname, buffer, (sizeof(fname) - strlen(SimDir)) );
	sprintf(buffer,".h5");
    strncat(fname, buffer, (sizeof(fname) - strlen(SimDir)) );
	
	printf("  filename= %s\n", fname);
	
	// hdf5 handles
    hid_t       file_id, dataset_psiRe, dataset_psiIm, dataset_a; // handles
    hid_t       dataspace_a, dataspace_psiRe, dataspace_psiIm; // handles
    hid_t       memspace; 
    herr_t      status;
    
    // memspace variables
    hsize_t  dimsm[4] = {N,N,N,2};        // memory space dimensions
    hsize_t  count_out[4]    = {N,N,N,1}; // size of the hyperslab in memory
    hsize_t  stride_out[4]   = {1,1,1,2}; // hyperslab stride in memory
    hsize_t  offset_outRe[4] = {0,0,0,0}; // hyperslab offset in memory
    hsize_t  offset_outIm[4] = {0,0,0,1}; // hyperslab offset in memory
	
	// dataset variables
	hsize_t adim[1] = {1};
	hsize_t dims[3] = {N, N, N};
	
	// create hdf5 file
	if (stat(fname, &st) == 0)
	  remove(fname);
	file_id = H5Fcreate(fname, H5F_ACC_TRUNC, H5P_DEFAULT, H5P_DEFAULT);


    // [0] save a
    dataspace_a = H5Screate_simple(1, adim, NULL);
	dataset_a = H5Dcreate(file_id, "/a", H5T_NATIVE_DOUBLE, dataspace_a,
			H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT);

    status = H5Dwrite(dataset_a, H5T_NATIVE_DOUBLE, H5S_ALL, H5S_ALL,
		      H5P_DEFAULT, &a);
	assert(status != -1);

	H5Dclose (dataset_a);
    H5Sclose (dataspace_a);

    
    // [1] save psiRe
    dataspace_psiRe = H5Screate_simple(3, dims, NULL);
    dataset_psiRe = H5Dcreate(file_id, "/psiRe", H5T_NATIVE_DOUBLE, dataspace_psiRe,
			H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT);
		
    memspace = H5Screate_simple(4,dimsm,NULL);   
    status = H5Sselect_hyperslab(memspace, H5S_SELECT_SET, offset_outRe, stride_out, 
				 count_out, NULL);
			
	status = H5Dwrite(dataset_psiRe, H5T_NATIVE_DOUBLE, memspace, H5S_ALL,
		      H5P_DEFAULT, psi[0]);
    
    H5Sclose(memspace);
    H5Dclose (dataset_psiRe);
    H5Sclose (dataspace_psiRe);
    
    // [2] save psiIm
    dataspace_psiIm = H5Screate_simple(3, dims, NULL);
    dataset_psiIm = H5Dcreate(file_id, "/psiIm", H5T_NATIVE_DOUBLE, dataspace_psiIm,
			H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT);
		
    memspace = H5Screate_simple(4,dimsm,NULL);   
    status = H5Sselect_hyperslab(memspace, H5S_SELECT_SET, offset_outIm, stride_out, 
				 count_out, NULL);
			
	status = H5Dwrite(dataset_psiIm, H5T_NATIVE_DOUBLE, memspace, H5S_ALL,
		      H5P_DEFAULT, psi[0]);
    
    H5Sclose(memspace);
    H5Dclose (dataset_psiIm);
    H5Sclose (dataspace_psiIm);

   // close file
    H5Fclose (file_id);

	printf( "snap saved.\n\n" );	
}


int set_snap(double a) {
	int snap;
	snap = round( log10(a) * ((double) NOUT / LOGAFINAL) );
	return snap;
}
