/*
******************************************************
This file is the multi-GPU version of 2D Heat Equation 
using CUDA model. This implementation is based on the 
CPU version from
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

#include<stdio.h>
#include<stdlib.h>
#include<string.h>
#include<sys/time.h>
#include <fcntl.h>
#include <unistd.h>
#include<cuda.h>

#define I2D(ni, i, j) ((i) + (ni)*(j))
#define THREADS 128

int NUM_GPUS;
int gangs[3];
int vectors[3];
char cu_filename[2][512];

int main(int argc, char** argv)
{
    if(argc < 6)
    {
        printf("Usage: %s <num GPUs> <ni> <nj> <nstep> <output file>\n", argv[0]);
        exit(1);        
    }

    int ni, nj, nstep;
    double tfac, *temp1_h, *temp2_h, *temp_tmp;
    int i, j, i2d, istep;
    double temp_bl, temp_br, temp_tl, temp_tr;
	struct timeval tim;
	double tick, tock;
    FILE *fp;
    FILE *fp2;
    int file_size;
    
	int fd;
    int d;
    int num_dev;
   
    // domain size and number of timesteps (iterations)
    NUM_GPUS = atoi(argv[1]);
    ni = atoi(argv[2]);
    nj = atoi(argv[3]);
    nstep = atoi(argv[4]);
    
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
	tick = tim.tv_sec + (tim.tv_usec/1000000.0);
    

    double **dev_temp1;
    double **dev_temp2;
    CUdevice device;
    CUcontext *contexts;
    CUmodule *cu_module;
    CUfunction *cu_function;
    CUstream *streams;
    int rows, LDA;
    double bc_size;
    int d1, d2;

    bc_size = (ni+2)*sizeof(double);
    num_dev = NUM_GPUS;
    dev_temp1 = (double**)malloc(num_dev*sizeof(double*));
    dev_temp2 = (double**)malloc(num_dev*sizeof(double*));

    contexts = (CUcontext*)malloc(num_dev*sizeof(CUcontext));
    streams = (CUstream*)malloc(num_dev*sizeof(CUstream));
    cu_module = (CUmodule*)malloc(num_dev*sizeof(CUmodule));
    cu_function= (CUfunction*)malloc(num_dev*sizeof(CUfunction));

    rows = nj/num_dev;
    LDA = ni + 2;


    cuInit(0);

    for(d=0; d<num_dev; d++)
    {
        cuDeviceGet(&device, d);
        cuCtxCreate(&contexts[d], 0, device);
    }

    cuMemHostRegister(temp1_h, (ni+2)*(nj+2)*sizeof(double), CU_MEMHOSTREGISTER_PORTABLE);
    cuMemHostRegister(temp2_h, (ni+2)*(nj+2)*sizeof(double), CU_MEMHOSTREGISTER_PORTABLE);
    
    for(d=0; d<num_dev; d++)
    {
        cuCtxPushCurrent(contexts[d]);
        cuMemAlloc((void**)&dev_temp1[d], (rows+2)*LDA*sizeof(double));
        cuMemAlloc((void**)&dev_temp2[d], (rows+2)*LDA*sizeof(double));
        cuStreamCreate(&streams[d], 0);
        cuCtxPopCurrent(&contexts[d]);
    }
    
    for(d=0; d<num_dev; d++)
    {
        cuCtxPushCurrent(contexts[d]);
        //cuMemcpyHtoDAsync(dev_temp1[d], temp1_h + d*rows*LDA, (rows+2)*LDA*sizeof(double), streams[d]);
        //cuMemcpyHtoDAsync(dev_temp2[d], temp2_h + d*rows*LDA, (rows+2)*LDA*sizeof(double), streams[d]);
        //cuMemcpyHtoD(dev_temp1[d], temp1_h + I2D(LDA, 0, d*rows), (rows+2)*LDA*sizeof(double));
        //cuMemcpyHtoD(dev_temp2[d], temp2_h + I2D(LDA, 0, d*rows), (rows+2)*LDA*sizeof(double));
        cuMemcpyHtoDAsync(dev_temp1[d], temp1_h + I2D(LDA, 0, d*rows), (rows+2)*LDA*sizeof(double), streams[d]);
        cuMemcpyHtoDAsync(dev_temp2[d], temp2_h + I2D(LDA, 0, d*rows), (rows+2)*LDA*sizeof(double), streams[d]);
        cuCtxPopCurrent(&contexts[d]);
    }

    for(istep = 0; istep < nstep; istep++)
    {
        for(d=0; d<num_dev; d++)
        {
            cuCtxPushCurrent(contexts[d]);
            d1 = ni+2;
            d2 = rows+2;

            if(strcmp(cu_filename[d], "heat_kernel.ptx") != 0)
            {
                char* ptx_source;
                strncpy(cu_filename[d], "heat_kernel.ptx", 512);
            
                fp = fopen(cu_filename[d], "rb");
                fseek(fp, 0, SEEK_END);
                file_size = ftell(fp);
                ptx_source = (char*)malloc((file_size+1)*sizeof(char));
                fseek(fp, 0, SEEK_SET);
                fread(ptx_source, sizeof(char), file_size, fp);
                fclose(fp);
                ptx_source[file_size] = '\0';
                cuModuleLoadData(&cu_module[d], ptx_source);
                free(ptx_source);
                cuModuleGetFunction(&cu_function[d], cu_module[d], "step_kernel");
             //   printf("if read kernel file from %d\n", d);
            }else
            {
                cuModuleGetFunction(&cu_function[d], cu_module[d], "step_kernel");
             //   print("else read kernel file from %d\n", d);
            }
        
            void *args[] = {&d1, &d2, &tfac, &dev_temp1[d], &dev_temp2[d]};
            vectors[0] = THREADS;  vectors[1] = 1; vectors[2] = 1;
            gangs[0] = (ni+THREADS-1)/THREADS;  gangs[1] = nj; gangs[2] = 1;

            cuLaunchKernel(cu_function[d], gangs[0], gangs[1], gangs[2],
                            vectors[0], vectors[1], vectors[2],
                            0,
                            streams[d], args, NULL);
            cuCtxPopCurrent(&contexts[d]);
        }

        
        for(d=0; d<num_dev; d++)
        {
            cuCtxPushCurrent(contexts[d]);
            cuStreamSynchronize(streams[d]);
            cuCtxPopCurrent(&contexts[d]);
        }
        
        for(d=0; d<num_dev; d++)
        {
            cuCtxPushCurrent(contexts[d]);
            if(d > 0)
                cuMemcpyPeerAsync(dev_temp2[d], contexts[d], dev_temp2[d-1] + I2D(LDA, 0, rows), contexts[d-1], bc_size, streams[d]);
                //cuMemcpyPeer(dev_temp2[d], contexts[d], dev_temp2[d-1] + I2D(LDA, 0, rows), contexts[d-1], bc_size);
            if(d < num_dev - 1)
                cuMemcpyPeerAsync(dev_temp2[d] + I2D(LDA, 0, rows+1), contexts[d], dev_temp2[d+1] + I2D(LDA, 0, 1), contexts[d+1], bc_size, streams[d]);
                //cuMemcpyPeer(dev_temp2[d] + I2D(LDA, 0, rows+1), contexts[d], dev_temp2[d+1] + I2D(LDA, 0, 1), contexts[d+1], bc_size);
            cuCtxPopCurrent(&contexts[d]);
        }
        
        for(d=0; d<num_dev; d++)
        {
            cuCtxPushCurrent(contexts[d]);
            cuStreamSynchronize(streams[d]);
            cuCtxPopCurrent(&contexts[d]);
        }
        
        for(d=0; d<num_dev; d++)
        {
            temp_tmp = dev_temp1[d];
            dev_temp1[d] = dev_temp2[d];
            dev_temp2[d] = temp_tmp;
        }
    }

    for(d=0; d<num_dev; d++)
    {
        cuCtxPushCurrent(contexts[d]);
        cuStreamSynchronize(streams[d]);
        cuMemcpyDtoHAsync(temp1_h + I2D(LDA, 0, d*rows+1), dev_temp1[d] + I2D(LDA, 0, 1), rows*LDA*sizeof(double), streams[d]);
        cuCtxPopCurrent(&contexts[d]);
    }

    
    for(d=0; d<num_dev; d++)
    {
        cuCtxPushCurrent(contexts[d]);
        cuStreamSynchronize(streams[d]);
        cuMemFree(dev_temp1[d]);
        cuMemFree(dev_temp2[d]);
        cuCtxPopCurrent(&contexts[d]);
    }


    free(dev_temp1);
    free(dev_temp2);
    
    cuMemHostUnregister(temp1_h); 
    cuMemHostUnregister(temp2_h); 
	
    gettimeofday(&tim, NULL);
	tock = tim.tv_sec + (tim.tv_usec/1000000.0);
    printf("Execution time is: %.2f s\n", tock - tick);

/*    
    fp2 = fopen(argv[5], "w");
    fprintf(fp2, "%d %d\n", ni, nj);
    for (j=0; j < nj; j++) {
        for (i=0; i < ni; i++) {
            fprintf(fp2, "%.4f\n", temp1_h[i + ni*j]);
        }
    }
    fclose(fp2);
  */
     fd = creat(argv[5], 00666);
     fd = open(argv[5], O_WRONLY);
     write(fd, temp1_h, (size_t)(ni+2)*(nj+2)*sizeof(double));
     close(fd);   
}
