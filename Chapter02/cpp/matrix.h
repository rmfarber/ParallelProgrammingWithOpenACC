#pragma once

#include<cstdlib>

struct matrix {
  unsigned int num_rows;
  unsigned int nnz;
  unsigned int *row_offsets;
  unsigned int *cols;
  double *coefs;
};


void allocate_3d_poission_matrix(matrix &A, int N) {
  int num_rows=(N+1)*(N+1)*(N+1);
  int nnz=27*num_rows;
  A.num_rows=num_rows;
  A.row_offsets=(unsigned int*)malloc((num_rows+1)*sizeof(unsigned int));
  A.cols=(unsigned int*)malloc(nnz*sizeof(unsigned int));
  A.coefs=(double*)malloc(nnz*sizeof(double));

  int offsets[27];
  double coefs[27];
  int zstride=N*N;
  int ystride=N;
  
  int i=0;
  for(int z=-1;z<=1;z++) {
    for(int y=-1;y<=1;y++) {
      for(int x=-1;x<=1;x++) {
        offsets[i]=zstride*z+ystride*y+x;
        if(x==0 && y==0 && z==0)
          coefs[i]=27;
        else
          coefs[i]=-1;
        i++;
      }
    }
  }

  nnz=0;
  for(int i=0;i<num_rows;i++) {
    A.row_offsets[i]=nnz;
    for(int j=0;j<27;j++) {
      int n=i+offsets[j];
      if(n>=0 && n<num_rows) {
        A.cols[nnz]=n;
        A.coefs[nnz]=coefs[j];
        nnz++;
      }
    }
  }

  A.row_offsets[num_rows]=nnz;
  A.nnz=nnz;
}

void free_matrix(matrix &A) {
  unsigned int *row_offsets=A.row_offsets;
  unsigned int * cols=A.cols;
  double * coefs=A.coefs;

  free(row_offsets);
  free(cols);
  free(coefs);
}


