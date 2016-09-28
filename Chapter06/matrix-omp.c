#include <sys/time.h>
#include <stdio.h>

/* matrix-omp.c */
#define SIZE 4000
float a[SIZE][SIZE];
float b[SIZE][SIZE];
float c[SIZE][SIZE];

int main()
{
  int i,j,k;
  struct timeval tim;
  double t1, t2;
  
  // Initialize matrices.
  for (i = 0; i < SIZE; ++i) {
    for (j = 0; j < SIZE; ++j) {
      a[i][j] = (float)i + j;
      b[i][j] = (float)i - j;
      c[i][j] = 0.0f;
    }
  }

  gettimeofday(&tim, NULL);
  t1=tim.tv_sec+(tim.tv_usec/1000000.0);

  // Compute matrix multiplication.
#pragma omp parallel for default(none) shared(a,b,c) private(i,j,k)
  for (i = 0; i < SIZE; ++i) {
    for (j = 0; j < SIZE; ++j) {
      for (k = 0; k < SIZE; ++k) {
	c[i][j] += a[i][k] * b[k][j];
      }
    }
  }
  gettimeofday(&tim, NULL);
  t2=tim.tv_sec+(tim.tv_usec/1000000.0);
  printf("%.6lf seconds elapsed\n", t2-t1);

  return 0;
}
