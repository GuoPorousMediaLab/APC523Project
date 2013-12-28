#include <iostream>
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
//		U1 = cell_interfaces_[0]->get_cell(0)->get_U();	// read value from left interface's left cell
//		U2 = U_;	// read value from current cell
//		U3 = cell_interfaces_[1]->get_cell(1)->get_U();	// read value from right interface's right cell
//		
//		for (i = 0; i < 4; i++)
//		{
//			deltal = U2[i] - U1[i];
//			deltar = U3[i] - U2[i];
//			if (deltal * deltar > 0)
//			{
//				if (deltal > 0)
//					delta = fmin(deltal, deltar);
//				else
//					delta = fmax(deltal, deltar);
//			}
//			else
//				delta = 0.0;
//				
//			Ul[i] = U2[i] - 0.5 * delta;
//			Ur[i] = U2[i] + 0.5 * delta;
//		}
//		cell_interfaces_[0]->set_U2(Ul);
//		cell_interfaces_[1]->set_U1(Ur);

//	set to piecewise constant for testing 1st order method
		cell_interfaces_[0]->set_U2(U_);
		cell_interfaces_[1]->set_U1(U_);
	}
}
