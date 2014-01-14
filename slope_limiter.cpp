#include <iostream>
#include <cmath>
#include "cell.h"
#include "interface.h"

void Cell::reconstruct(int slope_limiter, char direction)
// Piecewise linear reconstruction
// U1, U2, U_ stores the values in the left, middle, right cells; Ul, Ur stores the values after reconstruction
// Reference: Toro, section 14.4
{
	double U1[4], U2[4], Ul[4], Ur[4];
	int i;
	int interfaceid1, interfaceid2, cellid1, cellid2;

	if (direction == 'x')
	// if solving in x direction, get the left and right interfaces
	{
		interfaceid1 = 0;
		interfaceid2 = 1;
	}
	else if (direction == 'y')
	// if solving in y direction, fet the bottom and top interfaces
	{
		interfaceid1 = 2;
		interfaceid2 = 3;
	}
	
	// get the cellid of the two neighboring cells
	cellid1 = cell_interfaces_[interfaceid1]->get_cellid(0);
	cellid2 = cell_interfaces_[interfaceid2]->get_cellid(1);

//==============================set up U1[4]==============================
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

//==============================set up U2[4]==============================
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

	// if in y direction, swap u and v
	double temp;
	if (direction == 'y')
	{
		temp = U1[1];
		U1[1] = U1[2];
		U1[2] = temp;
		
		temp = U2[1];
		U2[1] = U2[2];
		U2[2] = temp;
		
		temp = U_[1];
		U_[1] = U_[2];
		U_[2] = temp;
	}
	
//==============================apply slope limiter==============================
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
	
// if in y direction, swap back u and v
	if (direction == 'y')
	{
		temp = Ul[1];
		Ul[1] = Ul[2];
		Ul[2] = temp;
		
		temp = Ur[1];
		Ur[1] = Ur[2];
		Ur[2] = temp;
		
		temp = U_[1];
		U_[1] = U_[2];
		U_[2] = temp;	
	}

//==============================assign reconstructed values	
	// set up the interface states for the Riemann solver
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
	double p, u, v, h, a;
	double deltal[4], deltar[4], delta[4], dU[4];
	int i;
	
	// parameters of the middle cell
	p = (GAMMA - 1.0) * (U2[3] - 0.5 * (U2[1] * U2[1] + U2[2] * U2[2]) / U2[0]);
	h = (U2[3] + p) / U2[0];
	u = U2[1] / U2[0];
	v = U2[2] / U2[0];
	a = sqrt(GAMMA * p / U2[0]);

	// decompose the two jumps to characterstics of the middle cell	
	deltal[2] = (U2[2] - U1[2]) - v * (U2[0] - U1[0]);
	deltal[1] = (GAMMA - 1.0) / (a * a) * ((U2[0] - U1[0]) * (h - u * u) + u * (U2[1] - U1[1]) - (U2[3] - U1[3]) + v * ((U2[2] - U1[2]) - v * (U2[0] - U1[0])));
	deltal[0] = ((u + a) * (U2[0] - U1[0]) - (U2[1] - U1[1]) - a * deltal[1]) / (2.0 * a);
	deltal[3] = (U2[0] - U1[0]) - deltal[0] - deltal[1];
	
	deltar[2] = (U3[2] - U2[2]) - v * (U3[0] - U2[0]);
	deltar[1] = (GAMMA - 1.0) / (a * a) * ((U3[0] - U2[0]) * (h - u * u) + u * (U3[1] - U2[1]) - (U3[3] - U2[3]) + v * ((U3[2] - U2[2]) - v * (U3[0] - U2[0])));
	deltar[0] = ((u + a) * (U3[0] - U2[0]) - (U3[1] - U2[1]) - a * deltar[1]) / (2.0 * a);
	deltar[3] = (U3[0] - U2[0]) - deltar[0] - deltar[1];
	
	// minbee function
	for (i = 0; i < 4; i++)
	{
		if (deltal[i] * deltar[i] > 0)	// if not a local extremum, set slope to minmod(left slope, right slope)
		{
			if (deltal[i] > 0)
				delta[i] = fmin(deltal[i], deltar[i]);
			else
				delta[i] = fmax(deltal[i], deltar[i]);
		}
		else	// if a local extremum, set slope to zero
			delta[i] = 0.0;
	}
	
	// combine the jumps in charasteristics back into a jump in U vector
	dU[0] = delta[0] + delta[1] + delta[3];
	dU[1] = delta[0] * (u - a) + delta[1] * u + delta[3] * (u + a);
	dU[2] = delta[0] * v + delta[1] * v + delta[2] + delta[3] * v;
	dU[3] = delta[0] * (h - u * a) + delta[1] * 0.5 * (u * u + v * v) + delta[2] * v + delta[3] * (h + u * a);
	
	// calculate the values at cell boundaries by linear extrapolation
	for (i = 0; i < 4; i++)
	{
		Ul[i] = U2[i] - 0.5 * dU[i];
		Ur[i] = U2[i] + 0.5 * dU[i];
	}
}

void Cell::superbee(double* U1, double* U2, double* U3, double* Ul, double* Ur)
// superbee slope limiter
{
	double p, u, v, h, a;
	double deltal[4], deltar[4], delta[4], dU[4];
	int i;
	
	// parameters of the middle cell
	p = (GAMMA - 1.0) * (U2[3] - 0.5 * (U2[1] * U2[1] + U2[2] * U2[2]) / U2[0]);
	h = (U2[3] + p) / U2[0];
	u = U2[1] / U2[0];
	v = U2[2] / U2[0];
	a = sqrt(GAMMA * p / U2[0]);
	
	// decompose the two jumps to characterstics of the middle cell
	deltal[2] = (U2[2] - U1[2]) - v * (U2[0] - U1[0]);
	deltal[1] = (GAMMA - 1.0) / (a * a) * ((U2[0] - U1[0]) * (h - u * u) + u * (U2[1] - U1[1]) - (U2[3] - U1[3]) + v * ((U2[2] - U1[2]) - v * (U2[0] - U1[0])));
	deltal[0] = ((u + a) * (U2[0] - U1[0]) - (U2[1] - U1[1]) - a * deltal[1]) / (2.0 * a);
	deltal[3] = (U2[0] - U1[0]) - deltal[0] - deltal[1];
	
	deltar[2] = (U3[2] - U2[2]) - v * (U3[0] - U2[0]);
	deltar[1] = (GAMMA - 1.0) / (a * a) * ((U3[0] - U2[0]) * (h - u * u) + u * (U3[1] - U2[1]) - (U3[3] - U2[3]) + v * ((U3[2] - U2[2]) - v * (U3[0] - U2[0])));
	deltar[0] = ((u + a) * (U3[0] - U2[0]) - (U3[1] - U2[1]) - a * deltar[1]) / (2.0 * a);
	deltar[3] = (U3[0] - U2[0]) - deltar[0] - deltar[1];
	
	// superbee function	
	for (int i = 0; i < 4; i++)
	{
		if (deltal[i] * deltar[i] > 0)
		// if not a local extremum,
		//set slope to maxmod(minmod(left slope, 2*right slope), minmod(2*left slope, right slope))
		{
			if (deltal[i] > 0)
				delta[i] = fmax(fmin(2.0 * deltal[i], deltar[i]), fmin(deltal[i], 2.0 * deltar[i]));
			else
				delta[i] = fmin(fmax(2.0 * deltal[i], deltar[i]), fmax(deltal[i], 2.0 * deltar[i]));
		}
		else	// if a local extremum, set slope to zero
			delta[i] = 0.0;
	}
	
	// combine the jumps in charasteristics back into a jump in U vector
	dU[0] = delta[0] + delta[1] + delta[3];
	dU[1] = delta[0] * (u - a) + delta[1] * u + delta[3] * (u + a);
	dU[2] = delta[0] * v + delta[1] * v + delta[2] + delta[3] * v;
	dU[3] = delta[0] * (h - u * a) + delta[1] * 0.5 * (u * u + v * v) + delta[2] * v + delta[3] * (h + u * a);
	
	// calculate the values at cell boundaries by linear extrapolation
	for (i = 0; i < 4; i++)
	{
		Ul[i] = U2[i] - 0.5 * dU[i];
		Ur[i] = U2[i] + 0.5 * dU[i];
	}
}
