#include <iostream>
#include <cmath>
#include "cell.h"

void Cell::minbeex()
// Minbee slope limiter
{
	double *U1, *U2, *U3;
	double Ul[4], Ur[4];	
	double deltal, deltar, delta;
	int i;
	
	if (cell_interfaces_[0]->get_cell(0) == NULL)	// left interface is left boundary
	{
		cell_interfaces_[0]->set_U1(U_);
		cell_interfaces_[0]->set_U2(U_);
		cell_interfaces_[1]->set_U1(U_);
	}
	else if (cell_interfaces_[1]->get_cell(1) == NULL)	// right interface is right boundary
	{
		cell_interfaces_[1]->set_U1(U_);
		cell_interfaces_[1]->set_U2(U_);
		cell_interfaces_[0]->set_U2(U_);
	}
	else
	{
		U1 = cell_interfaces_[0]->get_cell(0)->get_U();	// read value from left interface's left cell
		U2 = U_;	// read value from current cell
		U3 = cell_interfaces_[1]->get_cell(1)->get_U();	// read value from right interface's right cell
		
		for (i = 0; i < 4; i++)
		{
			deltal = U2[i] - U1[i];
			deltar = U3[i] - U2[i];
			if (deltal * deltar > 0)	// if not a local extremum, set slope to minmod(left slope, right slope)
			{
				if (deltal > 0)
					delta = fmin(deltal, deltar);
				else
					delta = fmax(deltal, deltar);
			}
			else	// if a local extremum, set slope to zero
				delta = 0.0;
				
			Ul[i] = U2[i] - 0.5 * delta;
			Ur[i] = U2[i] + 0.5 * delta;
		}
		
		// store the reconstructed value in the corresponding interfaces
		cell_interfaces_[0]->set_U2(Ul);
		cell_interfaces_[1]->set_U1(Ur);

//	set to piecewise constant for testing 1st order method
//		cell_interfaces_[0]->set_U2(U_);
//		cell_interfaces_[1]->set_U1(U_);
	}
}

void Cell::superbeex()
// superbee slope limiter
{
	double *U1, *U2, *U3;
	double Ul[4], Ur[4];	
	double deltal, deltar, delta;
	int i;
	
	if (cell_interfaces_[0]->get_cell(0) == NULL)	// left interface is left boundary
	{
		cell_interfaces_[0]->set_U1(U_);
		cell_interfaces_[0]->set_U2(U_);
		cell_interfaces_[1]->set_U1(U_);
	}
	else if (cell_interfaces_[1]->get_cell(1) == NULL)	// right interface is right boundary
	{
		cell_interfaces_[1]->set_U1(U_);
		cell_interfaces_[1]->set_U2(U_);
		cell_interfaces_[0]->set_U2(U_);
	}
	else
	{
		U1 = cell_interfaces_[0]->get_cell(0)->get_U();	// read value from left interface's left cell
		U2 = U_;	// read value from current cell
		U3 = cell_interfaces_[1]->get_cell(1)->get_U();	// read value from right interface's right cell
		
		for (i = 0; i < 4; i++)
		{
			deltal = U2[i] - U1[i];
			deltar = U3[i] - U2[i];
			if (deltal * deltar > 0)
			// if not a local extremum,
			//set slope to maxmod(minmod(left slope, 2*right slope), minmod(2*left slope, right slope))
			{
				if (deltal > 0)
					delta = fmax(fmin(2.0 * deltal, deltar), fmin(deltal, 2.0 * deltar));
				else
					delta = fmin(fmax(2.0 * deltal, deltar), fmax(deltal, 2.0 * deltar));
			}
			else	// if a local extremum, set slope to zero
				delta = 0.0;
				
			Ul[i] = U2[i] - 0.5 * delta;
			Ur[i] = U2[i] + 0.5 * delta;
		}
		
		// store the reconstructed value in the corresponding interfaces
		cell_interfaces_[0]->set_U2(Ul);
		cell_interfaces_[1]->set_U1(Ur);

//	set to piecewise constant for testing 1st order method
//		cell_interfaces_[0]->set_U2(U_);
//		cell_interfaces_[1]->set_U1(U_);
	}
}

void Cell::predictx(double dt)
{
	double *Ul, *Ur, Fl[4], Fr[4], pl, pr, hl, hr;
	int i;
	
	Ul = cell_interfaces_[0]->get_U2();
	Ur = cell_interfaces_[1]->get_U1();
	
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
	
	for (i = 0; i < 4; i++)
	{
		Ul[i] = Ul[i] + dt / dx_ * (Fl[i] - Fr[i]);
		Ur[i] = Ur[i] + dt / dx_ * (Fl[i] - Fr[i]);
	}
}
