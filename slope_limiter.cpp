#include <iostream>
#include <cmath>
#include "cell.h"
#include "interface.h"

void Cell::slopeLimiterx(int slope_limiter)
{
	double *U1, *U2, *U3;
	double ghost_U1[4], ghost_U2[4];
	double Ul[4], Ur[4];	
	int i;
	int temp_cellid1, temp_cellid2;

	temp_cellid1 = cell_interfaces_[0]->get_cellid(0);
	temp_cellid2 = cell_interfaces_[1]->get_cellid(1);

	if (temp_cellid1 < 0)	// left interface is left boundary
	{
		U2 = U_;	// read value from current cell
		U3 = cell_interfaces_[1]->get_cell(1)->get_U();	// read value from right interface's right cell

		if (temp_cellid1 == -1 || temp_cellid1 == -2)
		{
			cout << "temp_cellid1 " << temp_cellid1 << endl;
			cell_interfaces_[0]->set_U1(U_);
			cell_interfaces_[0]->set_U2(U_);
			cell_interfaces_[1]->set_U1(U_);
		}
		else if (temp_cellid1 == -3)
		{
			for (i = 0; i < 4; ++i)
			{
				if (i == 1)
				{
					ghost_U1[i] = - U2[i];
					ghost_U2[i] = - U3[i];
				}
				else { 
					ghost_U1[i] = U2[i];
					ghost_U2[i] = U3[i];
				}
			}
			switch (slope_limiter)
			{
				case 1:
					minbee(ghost_U1, ghost_U2, U2, Ul, Ur);
					cell_interfaces_[0]->set_U1(Ur);
					minbee(ghost_U2, U2, U3, Ul, Ur);
					cell_interfaces_[0]->set_U2(Ul);
					cell_interfaces_[1]->set_U1(Ur);
					break;
				case 2:
					superbee(ghost_U1, ghost_U2, U2, Ul, Ur);
					cell_interfaces_[0]->set_U1(Ur);
					superbee(ghost_U2, U2, U3, Ul, Ur);
					cell_interfaces_[0]->set_U2(Ul);
					cell_interfaces_[1]->set_U1(Ur);
					break;
				default:
					minbee(ghost_U1, ghost_U2, U2, Ul, Ur);
					cell_interfaces_[0]->set_U1(Ur);
					minbee(ghost_U2, U2, U3, Ul, Ur);
					cell_interfaces_[0]->set_U2(Ul);
					cell_interfaces_[1]->set_U1(Ur);
			}
		}
		
	}
	else if (temp_cellid2 < 0)	// right interface is right boundary
	{
		U1 = cell_interfaces_[0]->get_cell(0)->get_U();
		U2 = U_;	// read value from current cell
		U3 = cell_interfaces_[1]->get_cell(1)->get_U();	// read value from right interface's right cell

		if (temp_cellid2 == -1 || temp_cellid2 == -2)
		{
			cout << "temp_cellid2 " << temp_cellid2 << endl;
			cell_interfaces_[1]->set_U1(U_);
			cell_interfaces_[1]->set_U2(U_);
			cell_interfaces_[0]->set_U1(U_);
		}
		else if (temp_cellid2 == -3)
		{
			for (i = 0; i < 4; ++i)
			{
				if (i == 1)
				{
					ghost_U1[i] = - U2[i];
					ghost_U2[i] = - U1[i];
				}
				else { 
					ghost_U1[i] = U2[i];
					ghost_U2[i] = U1[i];
				}
			}
			switch (slope_limiter)
			{
				case 1:
					minbee(U2, ghost_U1, ghost_U2, Ul, Ur);
					cell_interfaces_[1]->set_U2(Ul);
					minbee(U1, U2, ghost_U1, Ul, Ur);
					cell_interfaces_[1]->set_U1(Ur);
					cell_interfaces_[0]->set_U2(Ul);
					break;
				case 2:
					superbee(U2, ghost_U1, ghost_U2, Ul, Ur);
					cell_interfaces_[1]->set_U2(Ul);
					superbee(U1, U2, ghost_U1, Ul, Ur);
					cell_interfaces_[1]->set_U1(Ur);
					cell_interfaces_[0]->set_U2(Ul);
					break;
				default:
					minbee(U2, ghost_U1, ghost_U2, Ul, Ur);
					cell_interfaces_[1]->set_U2(Ul);
					minbee(U1, U2, ghost_U1, Ul, Ur);
					cell_interfaces_[1]->set_U1(Ur);
					cell_interfaces_[0]->set_U2(Ul);
			}
		}		
	}
	else
	{
		U1 = cell_interfaces_[0]->get_cell(0)->get_U();	// read value from left interface's left cell
		U2 = U_;	// read value from current cell
		U3 = cell_interfaces_[1]->get_cell(1)->get_U();	// read value from right interface's right cell

		switch (slope_limiter)
		{
			case 1:
				minbee(U1, U2, U3, Ul, Ur);
				break;
			case 2:
				superbee(U1, U2, U3, Ul, Ur);
				break;
			default:
				minbee(U1, U2, U3, Ul, Ur);
		}
		cell_interfaces_[0]->set_U2(Ul);
		cell_interfaces_[1]->set_U1(Ur);
	}
}

void Cell::slopeLimitery(int slope_limiter)
{
	double *U1, *U2, *U3;
	double ghost_U1[4], ghost_U2[4];
	double Ub[4], Ut[4];	
	int i;
	int temp_cellid1, temp_cellid2;

	temp_cellid1 = cell_interfaces_[2]->get_cellid(0);
	temp_cellid2 = cell_interfaces_[3]->get_cellid(1);
	
	if (temp_cellid1 < 0)	// bottom interface is bottom boundary
	{
		U2 = U_;	// read value from current cell
		U3 = cell_interfaces_[3]->get_cell(1)->get_U();	// read value from top interface's above cell

		if (temp_cellid1 == -1 || temp_cellid1 == -2)
		{
			cell_interfaces_[2]->set_U1(U_);
			cell_interfaces_[2]->set_U2(U_);
			cell_interfaces_[3]->set_U1(U_);
		}
		else if (temp_cellid1 == -3)
		{
			for (i = 0; i < 4; ++i)
			{
				if (i == 2)
				{
					ghost_U1[i] = - U2[i];
					ghost_U2[i] = - U3[i];
				}
				else { 
					ghost_U1[i] = U2[i];
					ghost_U2[i] = U3[i];
				}
			}
			switch (slope_limiter)
			{
				case 1:
					minbee(ghost_U1, ghost_U2, U2, Ub, Ut);
					cell_interfaces_[2]->set_U1(Ut);
					minbee(ghost_U2, U2, U3, Ub, Ut);
					cell_interfaces_[2]->set_U2(Ub);
					cell_interfaces_[3]->set_U1(Ut);
					break;
				case 2:
					superbee(ghost_U1, ghost_U2, U2, Ub, Ut);
					cell_interfaces_[2]->set_U1(Ut);
					superbee(ghost_U2, U2, U3, Ub, Ut);
					cell_interfaces_[2]->set_U2(Ub);
					cell_interfaces_[3]->set_U1(Ut);
					break;
				default:
					minbee(ghost_U1, ghost_U2, U2, Ub, Ut);
					cell_interfaces_[2]->set_U1(Ut);
					minbee(ghost_U2, U2, U3, Ub, Ut);
					cell_interfaces_[2]->set_U2(Ub);
					cell_interfaces_[3]->set_U1(Ut);
			}
		}
		
	}
	else if (temp_cellid2 < 0)	// top interface is top boundary
	{
		U1 = cell_interfaces_[2]->get_cell(0)->get_U();

		if (temp_cellid2 == -1 || temp_cellid2 == -2)
		{
			cell_interfaces_[3]->set_U1(U_);
			cell_interfaces_[3]->set_U2(U_);
			cell_interfaces_[2]->set_U1(U_);
		}
		else if (temp_cellid2 == -3)
		{
			for (i = 0; i < 4; ++i)
			{
				if (i == 2)
				{
					ghost_U1[i] = - U2[i];
					ghost_U2[i] = - U1[i];
				}
				else { 
					ghost_U1[i] = U2[i];
					ghost_U2[i] = U1[i];
				}
			}
			switch (slope_limiter)
			{
				case 1:
					minbee(U2, ghost_U1, ghost_U2, Ub, Ut);
					cell_interfaces_[3]->set_U2(Ub);
					minbee(U1, U2, ghost_U1, Ub, Ut);
					cell_interfaces_[2]->set_U2(Ub);
					cell_interfaces_[3]->set_U1(Ut);
					break;
				case 2:
					superbee(U2, ghost_U1, ghost_U2, Ub, Ut);
					cell_interfaces_[3]->set_U2(Ub);
					superbee(U1, U2, ghost_U1, Ub, Ut);
					cell_interfaces_[2]->set_U2(Ub);
					cell_interfaces_[3]->set_U1(Ut);
					break;
				default:
					minbee(U2, ghost_U1, ghost_U2, Ub, Ut);
					cell_interfaces_[3]->set_U2(Ub);
					minbee(U1, U2, ghost_U1, Ub, Ut);
					cell_interfaces_[2]->set_U2(Ub);
					cell_interfaces_[3]->set_U1(Ut);
			}
		}		
	}
	else
	{
		U1 = cell_interfaces_[2]->get_cell(0)->get_U();	// read value from left interface's left cell
		U2 = U_;	// read value from current cell
		U3 = cell_interfaces_[3]->get_cell(1)->get_U();	// read value from right interface's right cell

		switch (slope_limiter)
		{
			case 1:
				minbee(U1, U2, U3, Ub, Ut);
				break;
			case 2:
				superbee(U1, U2, U3, Ub, Ut);
				break;
			default:
				minbee(U1, U2, U3, Ub, Ut);
		}
		cell_interfaces_[2]->set_U2(Ub);
		cell_interfaces_[3]->set_U1(Ut);
	}
}

void Cell::minbee(double* U1, double* U2, double* U3, double* Ul, double* Ur)
// Minbee slope limiter
{
	double deltal, deltar, delta;

	for (int i = 0; i < 4; i++)
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
}

void Cell::superbee(double* U1, double* U2, double* U3, double* Ul, double* Ur)
// superbee slope limiter
{
	double deltal, deltar, delta;	
	for (int i = 0; i < 4; i++)
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
}

void Cell::predictx(double dt)
{
	double *Ul, *Ur;
	Ul = cell_interfaces_[0]->get_U2();
	Ur = cell_interfaces_[1]->get_U1();	
	predict(Ul, Ur, dt);
}

void Cell::predicty(double dt)
{
	double *Ub, *Ut;
	Ub = cell_interfaces_[2]->get_U2();
	Ut = cell_interfaces_[3]->get_U1();
	predict(Ub, Ut, dt);
}

void Cell::predict(double* Ul, double* Ur, double dt)
{
	double Fl[4], Fr[4], pl, pr, hl, hr;
	int i;
	
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
