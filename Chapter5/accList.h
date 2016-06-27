#ifndef ACC_LIST_H
#define ACC_LIST_H

#include <cstdlib>
#include <cassert>
#ifdef _OPENACC
#include <openacc.h>
#endif

template<typename T>
class accList {
 public:
  explicit accList() {} 
  explicit accList(size_t size) {
    #pragma acc enter data copyin(this)
    allocate(size);
  }

  ~accList() {
    release();
    #pragma acc exit data delete(this)
  }

  #pragma acc routine seq
  T& operator[](size_t idx) { return _A[idx]; }

  #pragma acc routine seq
  const T& operator[](size_t idx) const { return _A[idx]; }

  size_t size() const { return _size; }

  accList& operator=(const accList& B) {
    allocate(B.size());
    for (size_t j = 0; j < _size; ++j) {
      _A[j] = B[j];
    }
    accUpdateDevice();
    return *this;
  }

  void insert(size_t idx, const T& val) { _A[idx] = val; }
  void insert(size_t idx, const T* val) { _A[idx] = *val; }

  void accUpdateSelf() { accUpdateSelfT(_A, 0); }
  void accUpdateDevice() { accUpdateDeviceT(_A, 0); }

 private:
  void release() {
    if (_size > 0) {
      #pragma acc exit data delete(_A[0:_size])
      delete[] _A;
      _A = nullptr;
      _size = 0;
    }
  }

  void allocate(size_t size) {
    if (_size != size) {
      release();
      _size = size;
      #pragma acc update device(_size)
      if (_size > 0) {
        _A = new T[_size];
#ifdef _OPENACC
        assert(!acc_is_present(&_A[0],sizeof(T)));
#endif
        #pragma acc enter data create(_A[0:_size])
      }
    }
  }

  // These template functions use the SFINAE pattern so as to call
  // the element type's accUpdateSelf/Device() member functions, if
  // they exist, on all of the elements, or otherwise use the
  // "acc update self/device" directive on the raw data.
  // This works because the versions of the template functions below 
  // that call those member functions look more specialized so 
  // are favored by C++'s function overload resolution rules  
  // over the more general versions that apply the directives.

  template<typename U>
  void accUpdateSelfT(U *p, long) {
    #pragma acc update self(p[0:_size])
  }

  template<typename U>
  auto accUpdateSelfT(U *p, int) -> decltype(p->accUpdateSelf()) {
    for (size_t j = 0; j < _size; ++j) {
      p[j].accUpdateSelf();
    }
  }

  template<typename U>
  void accUpdateDeviceT(U *p, long) {
    #pragma acc update device(p[0:_size])
  }

  template<typename U>
  auto accUpdateDeviceT(U *p, int) -> decltype(p->accUpdateDevice()) {
    for (size_t j = 0; j < _size; ++j) {
      p[j].accUpdateDevice();
    }
  }

  T* _A{nullptr};   // Contained objects
  size_t _size{0};  // Number of contained objects
};
#endif
