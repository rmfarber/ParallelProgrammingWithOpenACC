/*
******************************************************
This file is the OpenACC single GPU version of 2D Heat Equation 
using OpenACC programming model. This implementation is based on
the CPU version from
http://www.many-core.group.cam.ac.uk/archive/CUDAcourse09/

Permission to use, copy, distribute and modify this software for any 
purpose with or without fee is hereby granted. This software is        
provided "as is" without express or implied warranty. 

Send comments or suggestions for this OpenACC version to
            rxu6@uh.edu, schandra@udel.edu

Authors: Rengan Xu, Sunita Chandrasekaran

May 26th, 2016
******************************************************
*/

#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include <sys/time.h>
#include <fcntl.h>
#include <unistd.h>

// I2D to index into a linear memory space from a 2D array index pair
#define I2D(ni, i, j) ((i) + (ni)*(j))


// kernel to update temperatures - CPU version
void step_kernel_cpu(int ni, 
                     int nj,
                     double tfac, 
                     double *temp_in,
                     double *temp_out) {
    int i, j, i00, im10, ip10, i0m1, i0p1;
    double d2tdx2, d2tdy2;

    // loop over all points in domain (not boundary points)
#pragma acc kernels present(temp_in[0:ni*nj], temp_out[0:ni*nj]) 
{
	#pragma acc loop independent
    for (j=1; j < nj-1; j++) {
	#pragma acc loop independent
        for (i=1; i < ni-1; i++) {
	    // find indices into linear memory for central point and neighbors
            i00 = I2D(ni, i, j);
            im10 = I2D(ni, i-1, j);
            ip10 = I2D(ni, i+1, j);
            i0m1 = I2D(ni, i, j-1);
            i0p1 = I2D(ni, i, j+1);

	    // evaluate derivatives
            d2tdx2 = temp_in[im10] - 2*temp_in[i00] + temp_in[ip10];
            d2tdy2 = temp_in[i0m1] - 2*temp_in[i00] + temp_in[i0p1];

	    // update temperatures
            temp_out[i00] = temp_in[i00] + tfac*(d2tdx2 + d2tdy2);
        }
    }

}
}


int main(int argc, char *argv[]) 
{
    
    if(argc < 5)
    {
        printf("Usage: %s <ni> <nj> <nstep> <output file>\n", argv[0]);
        exit(0);
    }

    int ni, nj, nstep;
    double tfac, *temp1_h, *temp2_h, *temp_tmp;
    int i, j, i2d, istep;
    double temp_bl, temp_br, temp_tl, temp_tr;
	struct timeval tim;
	double start, end;
    double time;
	int fd;
   
    // domain size and number of timesteps (iterations)
    ni = atoi(argv[1]);
    nj = atoi(argv[2]);
    nstep = atoi(argv[3]);
    
    // allocate temperature array on host
    temp1_h = (double *)malloc(sizeof(double)*(ni+2)*(nj+2));
    temp2_h = (double *)malloc(sizeof(double)*(ni+2)*(nj+2));

    // initial temperature in interior
    for (j=1; j < nj+1; j++) {
        for (i=1; i < ni+1; i++) {
            i2d = i + (ni+2)*j;
            temp1_h[i2d] = 0.0;
        }
    }
    
    // initial temperature on boundaries - set corners
    temp_bl = 200.0f;
    temp_br = 300.0f;
    temp_tl = 200.0f;
    temp_tr = 300.0f;

    // set edges by linear interpolation from corners
    for (i=0; i < ni+2; i++) {
        // bottom
        j = 0;
        i2d = i + (ni+2)*j;
        temp1_h[i2d] = temp_bl + (temp_br-temp_bl)*(double)i/(double)(ni+1);

        // top
        j = nj+1;
        i2d = i + (ni+2)*j;
        temp1_h[i2d] = temp_tl + (temp_tr-temp_tl)*(double)i/(double)(ni+1);
    }

    for (j=0; j < nj+2; j++) {
        // left
        i = 0;
        i2d = i + (ni+2)*j;
        temp1_h[i2d] = temp_bl + (temp_tl-temp_bl)*(double)j/(double)(nj+1);

        // right
        i = ni+1;
        i2d = i + (ni+2)*j;
        temp1_h[i2d] = temp_br + (temp_tr-temp_br)*(double)j/(double)(nj+1);
    }

    // duplicate temeperature array on host
    memcpy(temp2_h, temp1_h, sizeof(double)*(ni+2)*(nj+2));
    
        
    tfac = 0.2;
    
	gettimeofday(&tim, NULL);
	start = tim.tv_sec + (tim.tv_usec/1000000.0);

    // main iteration loop  

#pragma acc data copy(temp1_h[0:(ni+2)*(nj+2)]) \
copyin(temp2_h[0:(ni+2)*(nj+2)]) deviceptr(temp_tmp)
{
    for (istep=0; istep < nstep; istep++) {
            // CPU kernel
            step_kernel_cpu(ni+2, nj+2, tfac, temp1_h, temp2_h);

	    // swap the temp pointers
            temp_tmp = temp1_h;
            temp1_h = temp2_h;
            temp2_h = temp_tmp;
            
    } 
}
	gettimeofday(&tim, NULL);
    end = tim.tv_sec + (tim.tv_usec/1000000.0);
    printf("Time for computing: %.2f s\n",end-start);


    // output temp1 to a file
    
	fd = creat(argv[4], 00666);
	fd = open(argv[4], O_WRONLY);
	write(fd, temp1_h, (size_t)(ni+2)*(nj+2)*sizeof(double));
	close(fd);
    
/*
    FILE *fp;
    fp = fopen(filename, "w");
    fprintf(fp, "%d %d\n", ni, nj);
    for (j=0; j < nj; j++) {
        for (i=0; i < ni; i++) {
            fprintf(fp, "%.4f\n", j, i, temp1_h[i + ni*j]);
        }
    }
    fclose(fp);
*/
}


