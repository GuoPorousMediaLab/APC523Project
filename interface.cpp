#include "interface.h"
#include <iostream>

void Interface::create_(vector<int> cell_id, int id, Model* mymodel)
{
	cellid_ = cell_id;
	id_ = id;
	mymodel_ = mymodel;
}

void Interface::create_()
{
	id_ = 0;
}

void Interface::set_interface_cells()
{
	vector<int>::iterator cellid_iterator = cellid_.begin();
	vector<int>::iterator cellid_end = cellid_.end();
	while (cellid_iterator != cellid_end)
	{
		interface_cells_.push_back(&(*mymodel_).get_cells()[*cellid_iterator - 1]);
		cellid_iterator++;
	}
}

void Interface::initialize()
{
	// to be implemented from IO 
}

void Interface::OutputCellid()
{
	vector<Cell*>::iterator interface_cells_iter = interface_cells_.begin();
    vector<Cell*>::iterator interface_cells_end = interface_cells_.end();
    while (interface_cells_iter != interface_cells_end) {
        cout << (*interface_cells_iter)->get_id() << endl;
        interface_cells_iter++;
    }
}