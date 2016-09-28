#include <stdio.h>
#include <math.h>
#include <stdlib.h>

#define NX 256
#define NSTEPS 500

int main(int argc, char **argv) {

   int i,n;

   int nx = NX;
   int nsteps = NSTEPS;
   int dummy;
   char refFile[64] = "wave3d.xline.ref";
   char newFile[64] = "wave3d.xline";

   if( argc >= 2 ) {
      nx = atoi( argv[1] );
   }
   if( argc >=3 ) {
      nsteps = atoi( argv[2] );
   }
   if( argc >= 5 ) {
      strcpy(refFile, argv[3]);
      strcpy(newFile, argv[4]);
   }

   size_t nbytes = nx * sizeof(float);

   printf("nx = %d, nsteps = %d, nbytes = %d\n", nx, nsteps, nbytes);

   float *restrict x = (float*)malloc( nbytes );
   float *restrict f = (float*)malloc( nbytes );
   float *restrict x_ref  = (float*)malloc( nbytes );
   float *restrict f_ref  = (float*)malloc( nbytes );
   if( 0==x || 0==f || 0==x_ref || 0==f_ref ) {
      printf( "couldn't allocate fields on the host\n" );
      return (-1);
   }

   float eps = 1.0e-08;
   float step = 0.0f;
   float dt = 0.000005f; // in order for the system to be numerically dt < dx!!!
   int OK = 1;
   int printevery = 20;
   
   FILE *fPtrRef = fopen(refFile, "r");
   for (i=0; i<nx; i++) {
      fscanf(fPtrRef,"%g %g",&x_ref[i], &f_ref[i]);
   }
   fscanf(fPtrRef,"%d", &dummy);
   
   FILE *fPtr = fopen(newFile, "r");
   for (i=0; i<nx; i++) {
      fscanf(fPtr,"%g %g",&x[i], &f[i]);
   }
   fscanf(fPtr,"%d", &dummy);
   
   for (i=0; i<nx; i++) {
      if (fabs(x_ref[i]-x[i]) >= eps || fabs(f_ref[i]-f[i]) >= eps) {
         printf("%5.3f, %5.3f, %10.6e, %10.6e\n", x_ref[i],x[i],f_ref[i],f[i]);
         OK = 0;
      }
   }
   
   for (n=0; n<nsteps; n++) {
      step = step + dt;
      if (((n+1)%printevery)==0) {
         printf("step = %9.6f \n",step);
         for (i=0; i<nx; i++) {
            fscanf(fPtrRef,"%g %g", &x_ref[i], &f_ref[i]);
         }
         fscanf(fPtrRef,"%d", &dummy);
   
         for (i=0; i<nx; i++) {
            fscanf(fPtr,"%g %g", &x[i], &f[i]);
         }
         fscanf(fPtr,"%d", &dummy);
         
         for (i=0; i<nx; i++) {
            if (fabs(x_ref[i]-x[i]) >= eps || fabs(f_ref[i]-f[i]) >= eps) {
               if ((fabs(x_ref[i]-x[i])/x_ref[i]) >= eps || fabs((f_ref[i]-f[i])/f_ref[i]) >= eps) {
                  printf("%d, %5.3f, %10.6e, %5.3f, %10.6e\n", n, x_ref[i], x[i], f_ref[i], f[i]);
                  OK = 0;
               }
            }     
         }
      }
   }
   
   if (OK) {
      printf("Results are identical\n");
   } else {
      printf("Results are not identical\n");
   }

   fclose(fPtrRef);
   fclose(fPtr);
   
   free(x);
   free(f);
   free(x_ref);
   free(f_ref);
   
   return(0);
}
