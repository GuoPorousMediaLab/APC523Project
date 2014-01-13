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

void Cell::initialize(double *U)
{
	for (int i = 0; i < 4; i++)
	{
		U_[i] = U[i];
	}
}

void Cell::predict(double dt, char direction)
{
	double Ul[4], Ur[4], Fl[4], Fr[4], pl, pr, hl, hr, temp;
	int interfaceid1, interfaceid2, i;

	if (direction == 'x')
	{
		interfaceid1 = 0;
		interfaceid2 = 1;
	}
	else if (direction == 'y')
	{
		interfaceid1 = 2;
		interfaceid2 = 3;
	}
	for (i = 0; i < 4; i++)	// get the reconstructed values Ul, Ur at cell boundaries
	{
		Ul[i] = cell_interfaces_[interfaceid1]->get_U2()[i];
		Ur[i] = cell_interfaces_[interfaceid2]->get_U1()[i];	
	}
	
	if (direction == 'y')
	{
		temp = Ul[1];
		Ul[1] = Ul[2];
		Ul[2] = temp;
		temp = Ur[1];
		Ur[1] = Ur[2];
		Ur[2] = temp;		
	}
	
	// calculate the flux Fl, Fr at the cell boundaries
	pl = (GAMMA - 1.0) * (Ul[3] - 0.5 * (Ul[1] * Ul[1] + Ul[2] * Ul[2]) / Ul[0]);
	pr = (GAMMA - 1.0) * (Ur[3] - 0.5 * (Ur[1] * Ur[1] + Ur[2] * Ur[2]) / Ur[0]);
	hl = (Ul[3] + pl) / Ul[0];
	hr = (Ur[3] + pr) / Ur[0];
	Fl[0] = Ul[1];
	Fl[1] = Ul[1] * Ul[1] / Ul[0] + pl;
	Fl[2] = Ul[1] * Ul[2] / Ul[0];
	Fl[3] = Ul[1] * hl;
	Fr[0] = Ur[1];
	Fr[1] = Ur[1] * Ur[1] / Ur[0] + pr;
	Fr[2] = Ur[1] * Ur[2] / Ur[0];
	Fr[3] = Ur[1] * hr;
	
	// evolve Ul and Ur
	if (direction == 'x')
	{
		for (i = 0; i < 4; i++)
		{
			Ul[i] = Ul[i] + dt / dx_ * (Fl[i] - Fr[i]);
			Ur[i] = Ur[i] + dt / dx_ * (Fl[i] - Fr[i]);
		}
	}
	else if (direction == 'y')
	{
		for (i = 0; i < 4; i++)
		{
			Ul[i] = Ul[i] + dt / dy_ * (Fl[i] - Fr[i]);
			Ur[i] = Ur[i] + dt / dy_ * (Fl[i] - Fr[i]);
		}
		temp = Ul[1];
		Ul[1] = Ul[2];
		Ul[2] = temp;
		temp = Ur[1];
		Ur[1] = Ur[2];
		Ur[2] = temp;
	}
	
	// set the evolved values to the interfaces
	cell_interfaces_[interfaceid1]->set_U2(Ul);
	cell_interfaces_[interfaceid2]->set_U1(Ur);
	
	// deal with reflective boundary
	if (cell_interfaces_[interfaceid1]->get_cellid(0) == -3)
	{
		if (direction == 'x')
		{
			Ul[1] = -Ul[1];
		}
		else if (direction == 'y')
		{
			Ul[2] = -Ul[2];
		}
		cell_interfaces_[interfaceid1]->set_U1(Ul);
	}
	if (cell_interfaces_[interfaceid2]->get_cellid(1) == -3)
	{
		if (direction == 'x')
		{
			Ur[1] = -Ur[1];
		}
		else if (direction == 'y')
		{
			Ur[2] = -Ur[2];
		}
		cell_interfaces_[interfaceid2]->set_U2(Ur);
	}
}

void Cell::update(double dt, char direction)
// update vector U using the flux F from the cell boundaries
{
	int i;
	
	if (direction == 'x')
	{
		for (i = 0; i < 4; i++)
		{
			U_[i] = U_[i] + dt / dx_  * (cell_interfaces_[0]->get_F()[i] - cell_interfaces_[1]->get_F()[i]);
		}
	}
	else if (direction == 'y')
	{
		for (i = 0; i < 4; i++)
		{
			U_[i] = U_[i] + dt / dy_  * (cell_interfaces_[2]->get_F()[i] - cell_interfaces_[3]->get_F()[i]);
		}
	}
}

double Cell::get_dt()
// calculate the maximum allowable step size (corresponding to CPL = 1)
{
	double p, u, a;
	
	u = sqrt(U_[1] * U_[1] + U_[2] * U_[2]) / U_[0];
	p = (GAMMA - 1.0) * (U_[3] - 0.5 * (U_[1] * U_[1] + U_[2] * U_[2]) / U_[0]);
	a = sqrt(GAMMA * p / U_[0]);
	
	return(fmin(dx_/ (u + a), dy_/ (u + a)));
}

void Cell::OutputInterfaceid()
{
    for (int i = 0; i < 4; i++)
    {
        cout << cell_interfaces_[i]->get_id() << endl;
    }
} 
