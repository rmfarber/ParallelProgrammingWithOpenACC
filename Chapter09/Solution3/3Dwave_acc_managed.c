#include <stdio.h>
#include <math.h>
#include <stdlib.h>
#include <cuda_profiler_api.h>
#include "timer.h"

#define NX 256
#define NY 256
#define NZ 256
#define NSTEPS 500

#define OFFSET(i, j, k, width, depth) ( (k) + (width) * ((j) + (i) * (depth)) )

//
// evolve a 3D scalar wave equation
// the two fields we'll evolve are:
//    f - the field
//    g - the time derivative of the field
//

int main(int argc, char *argv[]) {

   int i,j,k,n;

   int nx = NX;
   int ny = NY;
   int nz = NZ;
   int nsteps = NSTEPS;

   if( argc >= 4 ) {
      nx = atoi( argv[1] );
      ny = atoi( argv[2] );
      nz = atoi( argv[3] );
   }
   if( argc >=5 )
      nsteps = atoi( argv[4] );

   StartTimer();

   size_t nbytes = nx * ny * nz * sizeof(float);

   float *restrict x  = (float*)malloc( nbytes );
   float *restrict y  = (float*)malloc( nbytes );
   float *restrict z  = (float*)malloc( nbytes );
   float *restrict f  = (float*)malloc( nbytes );
   float *restrict g  = (float*)malloc( nbytes );
   float *restrict fp = (float*)malloc( nbytes );
   float *restrict gp = (float*)malloc( nbytes );
   if( 0==x || 0==y || 0==z || 0==f || 0==g || 0==fp || 0==gp ) {
      printf( "couldn't allocate fields on the host\n" );
      return (-1);
   }

   float dx = 2.0f/(nx-1);
   float dy = 2.0f/(ny-1);
   float dz = 2.0f/(nz-1);
   float dt = 0.00000005f; // in order for the system to be numerically dt < dx!!!

   // initialize the grid to run from -1 to 1 in each direction
   for (i=0; i<nx; i++) {
      for (j=0; j<ny; j++) {
         for (k=0; k<nz; k++) {
            int offset = OFFSET(i, j, k, ny, nz);
            x[offset] = -1.0f + (i)*dx;
            y[offset] = -1.0f + (j)*dy;
            z[offset] = -1.0f + (k)*dz;
         }
      }
   } 

   // initialize the field to be a gaussian
   for (i=0; i<nx; i++) {
      for (j=0; j<ny; j++) {
         for (k=0; k<nz; k++) {
            int offset = OFFSET(i, j, k, ny, nz);
            f[offset] = 0.2f * exp( - ( x[offset]*x[offset] + 
                                        y[offset]*y[offset] + 
                                        z[offset]*z[offset] ) / 0.05f);
            g[offset] = 0.0f;
         }
      }
   } 

   // output the initial data when there are an even number of points, 
   // pick a line closest to a coordinate axis
   FILE *fPtr = fopen("wave3d.xline", "w");
   for (i=0; i<nx; i++) {
      int offset = OFFSET(i, ny/2, nz/2, ny, nz);
      fprintf(fPtr,"%5.3f %10.6e\n",x[offset],f[offset]);
   }
   fprintf(fPtr,"\n");

   float step = 0.0f;
   int printevery = 20;
   printf("step = %9.6f \n",step);

   cudaProfilerStart();
   {

      for (n=0; n<nsteps; n++) {

         step = step + dt;
    
         if (((n+1)%printevery)==0)
            printf("step = %9.6f \n",step);
    
         #pragma acc kernels
         {

            // predictor
            #pragma acc loop independent collapse(2) gang
            for (i=0; i<nx; i++) {
               for (j=0; j<ny; j++) {
                  #pragma acc loop independent vector
                  for (k=0; k<nz; k++) {
                     int offset = OFFSET(i, j, k, ny, nz);
                     fp[offset] = f[offset] + dt*g[offset];
                  }
               }
            } 
      
            // static boundaries
            #pragma acc loop independent collapse(2)
            for (j=0; j<ny; j++) {
               for (k=0; k<nz; k++) {
                  int xbeg = OFFSET(0,    j, k, ny, nz);
                  int xend = OFFSET(nx-1, j, k, ny, nz);
                  gp[xbeg] = g[xbeg];
                  gp[xend] = g[xend];
               }
            } 
      
            #pragma acc loop independent collapse(2)
            for (i=0; i<nx; i++) {
               for (k=0; k<nz; k++) {
                  int ybeg = OFFSET(i,    0, k, ny, nz);
                  int yend = OFFSET(i, ny-1, k, ny, nz);
                  gp[ybeg] = g[ybeg];
                  gp[yend] = g[yend];
               }
            } 
      
            #pragma acc loop independent collapse(2)
            for (i=0; i<nx; i++) {
               for (j=0; j<ny; j++) {
                  int zbeg = OFFSET(i, j,    0, ny, nz); 
                  int zend = OFFSET(i, j, nz-1, ny, nz); 
                  gp[zbeg] = g[zbeg];
                  gp[zend] = g[zend];
               }
            } 
      
            // use the predictor to update gp
            #pragma acc loop independent collapse(2)
            for (i=1; i<nx-1; i++) {
               for (j=1; j<ny-1; j++) {
                  #pragma acc loop independent vector
                  for (k=1; k<nz-1; k++) {
                     int current = OFFSET(i, j, k, ny, nz);

                     int next_x = OFFSET(i+1,   j,   k,   ny, nz);
                     int next_y = OFFSET(i,   j+1,   k,   ny, nz);
                     int next_z = OFFSET(i,     j, k+1,   ny, nz);

                     int prev_x = OFFSET(i-1,   j,   k,   ny, nz);
                     int prev_y = OFFSET(i,   j-1,   k,   ny, nz);
                     int prev_z = OFFSET(i,     j, k-1,   ny, nz);

                     gp[current] = g[current] + dt * (
                                   (fp[next_x] - 2.0f * fp[current] + fp[prev_x]) / dx / dx +
                                   (fp[next_y] - 2.0f * fp[current] + fp[prev_y]) / dy / dy +
                                   (fp[next_z] - 2.0f * fp[current] + fp[prev_z]) / dz / dz );
                  }
               }
            } 
      
            // use the average g's to update f
            #pragma acc loop independent collapse(2)
            for (i=0; i<nx; i++) {
               for (j=0; j<ny; j++) {
                  #pragma acc loop independent vector
                  for (k=0; k<nz; k++) {
                     int offset = OFFSET(i, j, k, ny, nz);
                     fp[offset] = f[offset] + dt*(0.5f * (g[offset] + gp[offset]));
                  }
               }
            } 
      
            // now update all the variables
            #pragma acc loop independent collapse(2)
            for (i=0; i<nx; i++) {
               for (j=0; j<ny; j++) {
                  #pragma acc loop independent vector
                  for (int k=0; k<nz; k++) {
                     int offset = OFFSET(i, j, k, ny, nz);
                     f[offset] = fp[offset];
                     g[offset] = gp[offset];
                  }
               }
            } 
      
         } // pragma acc kernels
    
         if (((n+1)%printevery)==0) {
            #pragma acc update host(x[0:nx*ny*nz], f[0:nx*ny*nz])
            for (i=0; i<nx; i++) {
               int offset = OFFSET(i, ny/2, nz/2, ny, nz);
               fprintf(fPtr,"%5.3f %10.6e\n",x[offset],f[offset]);
            }
            fprintf(fPtr,"\n");
         }
    
    
      } // for nsteps

   } // pragma acc data

   cudaDeviceSynchronize();
   cudaProfilerStop();

   free(x);
   free(y);
   free(z);
   free(f);
   free(g);
   free(fp);
   free(gp);

   float totalTime = GetTimer();
   printf("Total time: %f seconds\n", totalTime / 1000.0f);

   exit(0);
}
