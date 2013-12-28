#include <iostream>
#include <cmath>
#include "cell.h"

using namespace std;

void Cell::create_(int *interfaceid, double x, double y, double dx, double dy, int id, Model *mymodel)
{
	for (int i = 0; i < 4; i++)
	{
		interfaceid_[i] = interfaceid[i];
	}
	x_ = x;
	y_ = y;
	dx_ = dx;
	dy_ = dy;
	id_ = id;
	mymodel_ = mymodel;
}

void Cell::create_()
{
	x_ = 0.0;
	y_ = 0.0;
	dx_ = 0.0;
	dy_ = 0.0;
}

void Cell::set_cell_interfaces()
{
	int i, id;
	for (i = 0; i < 4; i++)
	{
		id = interfaceid_[i];
		cell_interfaces_[i] = mymodel_->get_interface(id);
	}
}

void Cell::set_U(double *U)
{
	for (int i = 0; i < 4; i++)
	{
		U_[i] = U[i];
	}
}

void Cell::updatex(double dt)
// update vector U using the flux F of the left and right interfaces
{
	for (int i = 0; i < 4; i++)
	{
		U_[i] = U_[i] + dt / dx_  * (cell_interfaces_[0]->get_F()[i] - cell_interfaces_[1]->get_F()[i]);
	}
}

double Cell::get_dt()
// calculate the maximum allowable step size (corresponding to CPL = 1)
{
	double p, u, a;
	
	u = sqrt(U_[1] * U_[1] + U_[2] * U_[2]) / U_[0];
	p = (GAMMA - 1.0) * (U_[3] - 0.5 * (U_[1] * U_[1] + U_[2] * U_[2]) / U_[0]);
	a = sqrt(GAMMA * p / U_[0]);
	
	return(dx_/ (u + a));
}

void Cell::OutputInterfaceid()
{
    for (int i = 0; i < 4; i++)
    {
        cout << cell_interfaces_[i]->get_id() << endl;
    }
}    
