FC=pgfortran
FCFLAGS=-fast
LDFLAGS=${FCFLAGS}

cg.x: matrix.o vector.o main.o
	${FC} -o $@ ${LDFLAGS} $^

matrix.o: matrix.F90
vector.o: vector.F90
main.o: main.F90

.SUFFIXES: .o .F90
.F90.o:
	${FC} ${FCFLAGS} $< -c -o $@

.PHONY: clean
clean:
	rm -Rf cg.x pgprof* *.o *.mod core
