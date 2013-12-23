#include <vector>
#include <iostream>
#include "cell.h"

using namespace std;

void Cell::create_(vector<int> interfaceid, double x, double y, double dx, double dy, int id, Model* mymodel)
{
	interfaceid_ = interfaceid;
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
	vector<int>::iterator interfaceid_iterator = interfaceid_.begin();
	vector<int>::iterator interfaceid_end = interfaceid_.end();
	while (interfaceid_iterator != interfaceid_end)
	{
		cell_interfaces_.push_back(&(*mymodel_).get_interfaces()[*interfaceid_iterator - 1]);
		interfaceid_iterator++;
	}
}

void Cell::initialize()
{
	// to be implemented from IO
}

void Cell::OutputInterfaceid()
{
	vector<Interface*>::iterator cell_interfaces_iter = cell_interfaces_.begin();
    vector<Interface*>::iterator cell_interfaces_end = cell_interfaces_.end();
    while (cell_interfaces_iter != cell_interfaces_end) {
        cout << (*cell_interfaces_iter)->get_id() << endl;
        cell_interfaces_iter++;
    }  
}