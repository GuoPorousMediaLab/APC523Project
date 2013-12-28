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

void Interface::initialize()
{
	// to be implemented from IO 
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
