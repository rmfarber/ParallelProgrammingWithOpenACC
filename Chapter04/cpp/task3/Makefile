CXX=pgc++
CXXFLAGS=-fast -Minfo=all -ta=tesla

mandelbrot.x: mandelbrot.o main.o
	$(CXX) $(CXXFLAGS) -o $@ $^
mandelbrot_solution.x: mandelbrot.o main_solution.o
	$(CXX) $(CXXFLAGS) -o $@ $^

.PHONY: clean
clean:
	rm -f *.o core *.lst *.ptx *.pgprof *.cubin *.s *.x *.mod
