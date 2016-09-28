CXX=pgc++
CXXFLAGS=-fast -ta=tesla,lineinfo -Minfo=all -Mneginfo

cg.x: main.cpp matrix.h matrix_functions.h vector.h vector_functions.h
	${CXX} ${CXXFLAGS} main.cpp -o cg.x
clean:
	rm -Rf cg.x pgprof* core *.o
