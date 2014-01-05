OBJS = main.o model.o cell.o interface.o riemann_solver.o slope_limiter.o
CXX = g++  
CXXFLAGS = -O3 -g
LINKFLAGS = -O3 -g
Simulation: $(OBJS)
	$(CXX) $(OBJS) $(LINKFLAGS) -o $@
clean: 
	rm -f Simulation *.o *.txt *.csv *~
depend: 
	$(CXX) -MM $(CXXFLAGS) *.cpp > .depend

.cpp.o: $<
	$(CXX) $(CXXFLAGS) -c $<

-include .depend