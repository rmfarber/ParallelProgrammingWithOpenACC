
#ifdef _OPENACC
#include <openacc.h>
#endif

using namespace std;


template<typename T>
class myList {

   private:
      T* _A{nullptr};   
      size_t _size{0};  
  
   public:

    #pragma acc routine seq
    T& operator[](size_t idx) { return _A[idx]; };

    #pragma acc routine seq
    const T& operator[](size_t idx) const { return _A[idx]; };

    size_t size() const {
	return _size;
    }

    explicit myList() { }
    explicit myList(size_t size) {
	_size = size;
        _A = new T[_size];
	#pragma acc enter data copyin(this)
        #pragma acc enter data create(_A[0:_size])
    }

    ~myList() {
        #pragma acc exit data delete(_A[0:_size])
	#pragma acc exit data delete(this)
	delete [] _A;
	_A=NULL;
        _size=0;
    }

    inline void accUpdateSelf() {
        #pragma acc update self(_A[0:_size])
    } 
    inline void accUpdateDevice() {
        #pragma acc update device(_A[0:_size])
    } 
};    


