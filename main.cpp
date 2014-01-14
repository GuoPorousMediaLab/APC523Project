#include <iostream>
#include "model.h"

using namespace std;

int main()
{
//==============================User Defined Paramters==============================
	int	riemann_solver = 1,	// 1: Roe, 2: HLLE, 3: HLLC, default: Roe
		slope_limiter = 1;	// 1: MINBEE, 2: SUPERBEE, default: MINBEE
				
	double 	radius = 0.5,	// radius of the cyliner, typical value: 0.3-0.6
			mach = 3.0,	// Mach number, typical value: 3.0 - 10.0
			CPL = 0.9,	// CPL coefficient, typical value: 0.5-0.9
			tmax = 1.0,	// final time, typical value: 2.0-10.0
			presteps = 300;	// presmooth steps to generate initial condition
			
	int	division = 30;	// number of grids to approximate the cylinder, grid size = radius / division (in both x and y)
						// typical value: 20-60
	
	const char 	*gridfile = "grid.txt",
				*outputfile = "result.txt";
//==============================Other Parameters==============================
	double t, dt;
	int step;

//==============================Grid Generation==============================
	cout << "Generate grid. ";
	grid(radius, division, gridfile);
	cout << "Finished." << endl;
	
	cout << "Create cells and interfaces. ";
    Model mymodel(gridfile);
    cout << "Finished." << endl;

//==============================Initialization==============================  
	cout << "Initialize cells. ";
    mymodel.Initialize(mach);
    cout << "Finished." << endl;
	
	// Use the HLLE Riemann solver to presmooth the initial condition
	cout << "Evolve " << presteps << " steps to generate the initial condition." << endl;
	
    for (step = 0; step < presteps; step++)
    {
    	if (step % 50 == 0)
    	{
    		cout << "step = " << step << endl;
    	}
    	
    	// calculate dt using CPL condition
    	dt = mymodel.Timestep(CPL);
		
		// evolve dt/2 in x direction
    	mymodel.Reconstructx(1);
    	mymodel.Predictx(dt / 4.0);
    	mymodel.Riemannx(2);
    	mymodel.Updatex(dt / 2.0);
    	
    	// evolve dt in y direction
    	mymodel.Reconstructy(1);
    	mymodel.Predicty(dt / 2.0);
    	mymodel.Riemanny(2);
    	mymodel.Updatey(dt / 1.0);
    	
    	// evolve dt/2 in x direction
    	mymodel.Reconstructx(1);
    	mymodel.Predictx(dt / 4.0);
    	mymodel.Riemannx(2);
    	mymodel.Updatex(dt / 2.0);
    }
    cout << "Finished." << endl;

//==============================Calculation==============================
	// display some information
	cout << endl << "Start calculation." << endl;
	cout << "Riemann solver: ";
	switch (riemann_solver)
	{
		case 1:
			cout << "Roe; ";
			break;
		case 2:
			cout << "HLLE; ";
			break;
		case 3:
			cout << "HLLC; ";
			break;
		default:
			cout << "Roe; ";
	}
	cout << "Slope limiter: ";
	switch (slope_limiter)
	{
		case 1:
			cout << "MINBEE." << endl;
			break;
		case 2:
			cout << "SUPERBEE." << endl;
			break;
		default:
			cout << "MINBEE." << endl;
	}
	cout << "Termination time: " << tmax << " => Estimated steps: " << (int) (tmax / dt) << endl;
	
	
	// start calculation
	step = 0;
    for (t = 0; t < tmax; t = t + dt)
    {
    	if (step % 50 == 0)
    	{
    		cout << "step = " << step << ", t = " << t << endl;
    	}

    	// calculate dt using CPL condition    	
		dt = mymodel.Timestep(CPL);
		
		// evolve dt/2 in x direction		
    	mymodel.Reconstructx(slope_limiter);
    	mymodel.Predictx(dt / 4.0);
    	mymodel.Riemannx(riemann_solver);
    	mymodel.Updatex(dt / 2.0);

    	// evolve dt in y direction    	
    	mymodel.Reconstructy(slope_limiter);
    	mymodel.Predicty(dt / 2.0);
    	mymodel.Riemanny(riemann_solver);
    	mymodel.Updatey(dt / 1.0);

		// evolve dt/2 in x direction    	
    	mymodel.Reconstructx(slope_limiter);
    	mymodel.Predictx(dt / 4.0);
    	mymodel.Riemannx(riemann_solver);
    	mymodel.Updatex(dt / 2.0);
    	    	
    	step++;
    }
	cout << "step = " << step << ", t = " << t << endl;
	
//==============================Output==============================
    mymodel.Output(outputfile);
    cout << "Calculation completed. Results saved in " << outputfile << '.' << endl;
	
	return 0;
}
