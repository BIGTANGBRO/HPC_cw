CC       = g++
CXX      = g++
CXXFLAGS = -Wall -O2 -fopenmp
LDLIBS   = -lblas -lboost_program_options
HDRS = ReactionDiffusion.hpp
OBJS = main.o ReactionDiffusion.o

.PHONY: clean run

%.o : %.cpp $(HDRS)
	$(CXX) $(CXXFLAGS) -o $@ -c $<
	
main: $(OBJS)
	$(CXX) -fopenmp -o $@ $^ $(LDLIBS)
	
run: main
	./main

test1: main
	./main --Nx 101 --Ny 101 --a 0.75 --b 0.06 --eps 50.0 --mu1 5.0 --mu2 0.0 --c 1

test2: main
	./main --Nx 251 --Ny 251 --a 0.75 --b 0.06 --eps 13.0 --mu1 5.0 --mu2 0.0 --c 1

test3: main
	./main --Nx 101 --Ny 101 --a 0.5 --b 0.1 --eps 50.0 --mu1 5.0 --mu2 0.0 --c 1

test4: main
	./main --Nx 151 --Ny 81 --a 0.75 --b 0.0001 --eps 12.5 --mu1 1.0 --mu2 0.01 --c 1

clean:
	rm -f main *.o
