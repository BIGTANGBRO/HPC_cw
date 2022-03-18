CC       = g++
CXX      = g++
CXXFLAGS = -Wall -O2 -fopenmp
LDLIBS   = -lblas
HDRS = ReactionDiffusion.hpp
OBJS = main.o ReactionDiffusion.o

.PHONY: clean run

%.o : %.cpp $(HDRS)
	$(CXX) $(CXXFLAGS) -o $@ -c $<
	
main: $(OBJS)
	$(CXX) -fopenmp -o $@ $^ $(LDLIBS)
	
run: main
	./main

clean:
	rm -f main *.o
