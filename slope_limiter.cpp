#include <iostream>
#include <cmath>
#include "cell.h"
#include "interface.h"

void Cell::reconstruct(int slope_limiter, char direction)
{
	double U1[4], U2[4], Ul[4], Ur[4];
	// U1, U2 stores the values in the neighboring cells; Ul, Ur stores the values after reconstruction
	int i;
	int interfaceid1, interfaceid2, cellid1, cellid2;

	if (direction == 'x')
	// if solving in x direction, set the two interfaces as left and right
	{
		interfaceid1 = 0;
		interfaceid2 = 1;
	}
	else if (direction == 'y')
	// if solving in y direction, set the two interfaces as bottom and top
	{
		interfaceid1 = 2;
		interfaceid2 = 3;
	}
	
	// get the cellid of the two neighboring cells
	cellid1 = cell_interfaces_[interfaceid1]->get_cellid(0);
	cellid2 = cell_interfaces_[interfaceid2]->get_cellid(1);

	// set up U1[4]
	if (cellid1 > 0)
	// if left interface is not a boundary, copy left cell to U1
	{
		for (i = 0; i < 4; i++)
		{
			U1[i] = cell_interfaces_[interfaceid1]->get_cell(0)->get_U()[i];
		}
	}
	else if (cellid1 == -1)
	// if left interface is a fixed-flux boundary, create a ghost cell corresponding to that flux
	{
		for (i = 0; i < 4; i++)
		{
			U1[i] = cell_interfaces_[interfaceid1]->get_U1()[i];
		}
	}
	else
	// if left interface is a transmissive boundary, create a ghost cell identical to current cell
	{
		for (i = 0; i < 4; i++)
		{
			U1[i] = U_[i];
		}
		if (cellid1 == -3)
	// if left interface is a reflective boundary, create a ghost cell with (rho, -rho*u, rho*v, E) or (rho, rho*u, -rho*v, E)
		{
			if (direction == 'x')
			{
				U1[1] = -U_[1];
			}
			else if (direction == 'y');
			{
				U1[2] = -U_[2];
			}
		}
	}

	
	// set up U2[4]
	if (cellid2 > 0)
	// if right interface is not a boundary, copy right cell to U2
	{
		for (i = 0; i < 4; i++)
		{
			U2[i] = cell_interfaces_[interfaceid2]->get_cell(1)->get_U()[i];
		}
	}
	else if (cellid2 == -1)
	// if left interface is a fixed-flux boundary, create a ghost cell corresponding to that flux
	{
		for (i = 0; i < 4; i++)
		{
			U2[i] = cell_interfaces_[interfaceid2]->get_U2()[i];
		}
	}
	else
	// if right interface is a transmissive boundary, create a ghost cell identical to current cell
	{
		for (i = 0; i < 4; i++)
		{
			U2[i] = U_[i];
		}
		if (cellid2 == -3)
	// if right interface is a reflective boundary, create a ghost cell with (rho, -rho*u, rho*v, E) or (rho, rho*u, -rho*v, E)
		{
			if (direction == 'x')
			{
				U2[1] = -U_[1];
			}
			else if (direction == 'y')
			{
				U2[2] = -U_[2];
			}
		}
	}
	
	// apply slope limiter
	switch (slope_limiter)
	{
		case 1:
			minbee(U1, U_, U2, Ul, Ur);
			break;
		case 2:
			superbee(U1, U_, U2, Ul, Ur);
			break;
		default:
			minbee(U1, U_, U2, Ul, Ur);
	}
	
	// use the reconstructed value Ul and Ur to set up the interface states for the Riemann solver
	cell_interfaces_[interfaceid1]->set_U2(Ul);
	cell_interfaces_[interfaceid2]->set_U1(Ur);
	
	// set up U1 for left boundary
	if (cellid1 == -2)
	{
		cell_interfaces_[interfaceid1]->set_U1(Ul);
	}
	else if (cellid1 == -3)
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
	// set up U2 for right boundary
	if (cellid2 == -2)
	{
		cell_interfaces_[interfaceid2]->set_U2(Ur);
	}
	else if (cellid2 == -3)
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

void Cell::minbee(double* U1, double* U2, double* U3, double* Ul, double* Ur)
// minbee slope limiter
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
