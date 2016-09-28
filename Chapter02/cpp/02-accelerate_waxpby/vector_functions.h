/*
 *  Copyright 2014 NVIDIA Corporation
 *
 *  Licensed under the Apache License, Version 2.0 (the "License");
 *  you may not use this file except in compliance with the License.
 *  You may obtain a copy of the License at
 *
 *      http://www.apache.org/licenses/LICENSE-2.0
 *
 *  Unless required by applicable law or agreed to in writing, software
 *  distributed under the License is distributed on an "AS IS" BASIS,
 *  WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
 *  See the License for the specific language governing permissions and
 *  limitations under the License.
 */
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
  double *restrict xcoefs=x.coefs;
  double *ycoefs=y.coefs;
  double *wcoefs=w.coefs;

#pragma acc kernels copy(wcoefs[:w.n],ycoefs[0:y.n]) copyin(xcoefs[0:x.n])
  {
#pragma acc loop independent
    for(int i=0;i<n;i++) {
      wcoefs[i]=alpha*xcoefs[i]+beta*ycoefs[i];
    }
  }
}

