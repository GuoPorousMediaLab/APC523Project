#include <iostream>
#include <cstdlib>
#include "cell.h"
#include "interface.h"
#include "model.h"
#include "riemann_solver.h"
#include "slope_limiter.h"

using namespace std;

int main(int argc, char const *argv[])
{
	// 1 - Roe, 2 - HLLE, 3 - HLLC
	// 1 - minbee, 2 - superbee
	if (argc != 3)
	{
		cerr << "You must choose Riemann solver and slope-limiter, e.g.: " << endl;
		cerr << "./Simulation 1 1" << endl;
		abort();
	}

	int riemannSolver = atoi(argv[1]);
	int slopeLimiter = atoi(argv[2]);

	int Nx = 100;
    int Ny = 1;
    double t, dt, tmax = 0.2, CPL = 0.9;
    
    Model mymodel(Nx, Ny);
    
    mymodel.Initialize();
    
    int step = 1;
    for (t = 0; t < tmax; t = t + dt)
    {
    	cout << "step " << step << endl;
		dt = mymodel.Timestep(CPL);
    	mymodel.Reconstructx(slopeLimiter);
    	mymodel.Predictx(dt / 2.0);
    	mymodel.Riemannx(riemannSolver);
    	mymodel.Updatex(dt);
    	step++;
    }

    // for test
//    mymodel.Outputid();
    mymodel.Outputvalue("result.txt");
	
	return 0;
}
