#include "interface.h"
#include <iostream>
#include <cmath>

void Interface::create_(int *cellid, int id, Model *mymodel)
{
	for (int i = 0; i < 2; i++)
	{
		cellid_[i] = cellid[i];
	}
	id_ = id;
	mymodel_ = mymodel;
}

void Interface::create_()
{
	id_ = 0;
}

void Interface::set_interface_cells()
{
	int i, id;
	for (i = 0; i < 2; i++)
	{
		id = cellid_[i];
		if (id == -1)
		{
			interface_cells_[i] = NULL;
		}
		else
		{
			interface_cells_[i]= mymodel_->get_cell(id);
		}
	}
}

void Interface::initialize(char direction)
{
	double p, h, *U;
	int i;
	
	if (cellid_[0] != -1 && cellid_[1] != -1)	// if not a fixed-flux boundary, no need to initialize F_[4]
		return;
		
	if (cellid_[0] == -1)
	{
		U = interface_cells_[1]->get_U();
		for (i = 0; i < 4; i++)
		{
			U1_[i] = U[i];
		}
	}
	else
	{
		U = interface_cells_[0]->get_U();
		for (i = 0; i < 4; i++)
		{
			U2_[i] = U[i];
		}
	}
	
	p = (GAMMA - 1.0) * (U[3] - 0.5 * (U[1] * U[1] + U[2] * U[2]) / U[0]);
	h = (U[3] + p) / U[0];
	
	if (direction == 'x')
	{
		F_[0] = U[1];
		F_[1] = U[1] * U[1] / U[0] + p;
		F_[2] = U[1] * U[2] / U[0];
		F_[3] = U[1] * h;
	}
	else if (direction == 'y')
	{
		F_[0] = U[2];
		F_[1] = U[1] * U[2] / U[0];
		F_[2] = U[2] * U[2] / U[0] + p;
		F_[3] = U[2] * h;
	}	
}

void Interface::set_U1(double *U)
{
	for (int i = 0; i < 4; i++)
	{
		U1_[i] = U[i];
	}
}

void Interface::set_U2(double *U)
{
	for (int i = 0; i < 4; i++)
	{
		U2_[i] = U[i];
	}
}

void Interface::OutputCellid()
{
    for (int i = 0; i < 2; i++)
    {
    	if (interface_cells_[i] == NULL)
    	{
    		cout << -1 << endl;
    	}
    	else
    	{
			cout << interface_cells_[i]->get_id() << endl;
		}
    }
}
