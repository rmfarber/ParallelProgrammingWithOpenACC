
#ifndef FLOATS3_H
#define FLOATS3_H

#include <cstring>
#ifdef _OPENACC
#include <openacc.h>
#endif

using namespace std;

struct soaFloat3 {

  float *x{nullptr}, *y{nullptr}, *z{nullptr};

  soaFloat3() {
    // Don't "acc enter data copyin(this)" for the default 
    // constructor, so that an instantiation as part of 
    // an array, as in accList, retains address contiguity 
    // in both host and device address spaces.
  }

  explicit soaFloat3(size_t n) : _created{true} {
    #pragma acc enter data copyin(this)
    allocate(n);
  }

  ~soaFloat3() {
    release();
    if (_created) {
      #pragma acc exit data delete(this)
    }
  }

  soaFloat3& operator=(const soaFloat3 &B) {
    size_t Bsize{B.size()};
    if (Bsize != _size)
      allocate(Bsize);
    size_t bytes = _size * sizeof *x;
    memcpy(x, B.x, bytes);
    memcpy(y, B.y, bytes);
    memcpy(z, B.z, bytes);
#if defined(USE_ACC_MEMCPY) && defined(_OPENACC)
    acc_memcpy_device(acc_deviceptr(x), acc_deviceptr(B.x), bytes);
    acc_memcpy_device(acc_deviceptr(y), acc_deviceptr(B.y), bytes);
    acc_memcpy_device(acc_deviceptr(z), acc_deviceptr(B.z), bytes);
#else
    accUpdateDevice();
#endif
    return *this;
  }

#pragma acc routine vector
  void add(const soaFloat3& B) {
    #pragma acc loop vector
    for (size_t j = 0; j < _size; ++j) {
      x[j] += B.x[j];
      y[j] += B.y[j];
      z[j] += B.z[j];
    }
  }

#pragma acc routine seq
  void setValue(float xval, float yval, float zval) {
#pragma acc loop seq
    for (size_t j = 0; j < _size; ++j) {
      x[j] = xval;
      y[j] = yval;
      z[j] = zval;
    }
  }

  size_t size() const { return _size; }

  void accUpdateSelf() {
    #pragma acc update self(x[0:_size])
    #pragma acc update self(y[0:_size])
    #pragma acc update self(z[0:_size])
  }

  void accUpdateDevice() {
    #pragma acc update device(x[0:_size])
    #pragma acc update device(y[0:_size])
    #pragma acc update device(z[0:_size])
  }

 private:
  void release() {
    if (_size > 0) {
      #pragma acc exit data delete(x[0:_size])
      #pragma acc exit data delete(y[0:_size])
      #pragma acc exit data delete(z[0:_size])
      delete[] x;
      delete[] y;
      delete[] z;
      _size = 0;
      x = y = z = nullptr;
    }
  }

  void allocate(size_t n) {
    if (_size != n) {
      release();
      _size = n;
      #pragma acc update device(_size)
      if (_size > 0) {
        x = new float[_size];
        y = new float[_size];
        z = new float[_size];
        #pragma acc enter data create(x[0:_size])
        #pragma acc enter data create(y[0:_size])
        #pragma acc enter data create(z[0:_size])
      }
    }
  }

  size_t _size{0};
// Whether to "acc exit data delete(this)" in d'tor
  const bool _created{false};  
};
#endif
