#include <iostream>
#include "cell.h"
#include "interface.h"
#include "model.h"
#include "reiman_solver.h"
#include "slope_limiter.h"

using namespace std;

int main(int argc, char const *argv[])
{
	// 1 - Roe, 2 - HLLE, 3 - HLLC
	// 1 - minbee, 2 - superbee
	if (argc != 3)
	{
		cerr << "You must choose Rieman solver and slope limiter, e.g.: " << endl;
		cerr << "./Simulation 1 1" << endl;
		abort();
	}

	int riemanSolver = atoi(argv[1]);
	int slopeLimiter = atoi(argv[2]);

	int Nx = 4;
    int Ny = 4;
    Model mymodel(Nx, Ny);
    // for test
    mymodel.Outputid();
	
	return 0;
}