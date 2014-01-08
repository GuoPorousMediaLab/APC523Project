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

	int riemann_solver = atoi(argv[1]);
	int slope_limiter = atoi(argv[2]);

	int Nx = 100;
    int Ny = 1;
    double t, dt, tmax = 10.0, CPL = 0.9;
    
    Model mymodel(Nx, Ny);
    
    mymodel.Initialize();
    
    int step = 1;
    for (t = 0; t < tmax; t = t + dt)
    {
    	cout << "step " << step << endl;
		dt = mymodel.Timestep(CPL);
		
    	mymodel.Reconstructx(slope_limiter);
    	mymodel.Predictx(dt / 4.0);
    	mymodel.Riemannx(riemann_solver);
    	mymodel.Updatex(dt / 2.0);
    	
    	mymodel.Reconstructy(slope_limiter);
    	mymodel.Predicty(dt / 2.0);
    	mymodel.Riemanny(riemann_solver);
    	mymodel.Updatey(dt / 1.0);
    	
    	mymodel.Reconstructx(slope_limiter);
    	mymodel.Predictx(dt / 4.0);
    	mymodel.Riemannx(riemann_solver);
    	mymodel.Updatex(dt / 2.0);
    	
    	step++;
    }

    // for test
    //mymodel.Outputid();
    mymodel.Outputvalue("result.txt");
	
	return 0;
}
