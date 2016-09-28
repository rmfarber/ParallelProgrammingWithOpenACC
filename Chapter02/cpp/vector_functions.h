#pragma once
#include<cstdlib>
#include "vector.h"


double dot(const vector& x, const vector& y) {
  double sum=0;
  unsigned int n=x.n;
  double *xcoefs=x.coefs;
  double *ycoefs=y.coefs;

  for(int i=0;i<n;i++) {
    sum+=xcoefs[i]*ycoefs[i];
  }
  return sum;
}

void waxpby(double alpha, const vector &x, double beta, const vector &y, const vector& w) {
  unsigned int n=x.n;
  double *xcoefs=x.coefs;
  double *ycoefs=y.coefs;
  double *wcoefs=w.coefs;

  for(int i=0;i<n;i++) {
    wcoefs[i]=alpha*xcoefs[i]+beta*ycoefs[i];
  }
}

