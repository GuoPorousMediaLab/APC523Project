#include <vector>
#include <iostream>
#include "model.h"

using namespace std;

void Model::create_(int Nx, int Ny)
{
	Nx_ = Nx;
	Ny_ = Ny;

	// initialize cells
	vector<int> interfaceid; // interfaceid for each cell
	double centerx = 0.0;
    double centery = 0.0;
    double dx = 0.4;
    double dy = 0.4;
    double startx = centerx - (Nx_-1) / 2.0 * dx;
    double starty = centery - (Ny_-1) / 2.0 * dy;
    for (int i = 1; i <= Ny_; i++)
    {    
    	for (int j = 1; j <= Nx_; j++)
        {
        	interfaceid.push_back((i-1) * Ny_ + j);
        	interfaceid.push_back(i * Ny_ + j);
        	interfaceid.push_back(Nx_ * (Ny_ + 1) + (i - 1) * (Ny_ + 1) + j);
        	interfaceid.push_back(Nx_ * (Ny_ + 1) + (i - 1) * (Ny_ + 1) + j + 1);
            cells_.push_back(Cell(interfaceid, startx + (j-1)*dx, starty + (i-1)*dy, dx, dy, (i - 1) * Nx_ + j, this));
            interfaceid.clear();
        }
    }

	// initialize intefaces
	vector<int> cellid; // cellid for each interface
	for (int i = 1; i <= Ny_ + 1; ++i)
	{
		for (int j = 1; j <= Nx_; ++j)
		{
			if (i < Ny_ + 1) {
                cellid.push_back((i-1) * Nx_ + j);
            }
            if (i > 1) {
                cellid.push_back((i-2) * Nx_ + j);
            }
			interfaces_.push_back(Interface(cellid, (i - 1) * Nx_ + j, this));
            cellid.clear();
		}
	}

    for (int i = 1; i <= Ny_; ++i)
    {
        for (int j = 1; j <= Nx_ + 1; ++j)
        {
            if (j < Nx_ + 1) {
                cellid.push_back((i - 1) * Nx_ + j);
            }
            if (j > 1) {
                cellid.push_back((i - 1) * Nx_ + j - 1);
            }
            interfaces_.push_back(Interface(cellid, Nx_ * (Ny_ + 1) + (i - 1) * (Nx_ + 1) + j, this));
            cellid.clear();
        }
    }

    // setup cells correspond to interface and interfaces correspond to cells
    vector<Interface>::iterator interfaces_iter = interfaces_.begin();
    vector<Interface>::iterator interfaces_end = interfaces_.end();
    while (interfaces_iter != interfaces_end)
    {
        (*interfaces_iter).set_interface_cells();
        interfaces_iter++;
    }

    vector<Cell>::iterator cells_iter = cells_.begin();
    vector<Cell>::iterator cells_end = cells_.end();
    while (cells_iter != cells_end)
    {
        (*cells_iter).set_cell_interfaces();
        cells_iter++;
    }

}

void Model::create_()
{
	Nx_ = 0;
	Ny_ = 0;
}

void Model::Outputid()
{
    vector<Interface>::iterator interfaces_iter = interfaces_.begin();
    vector<Interface>::iterator interfaces_end  = interfaces_.end();
    while (interfaces_iter != interfaces_end)
    {
        cout << "interfaceid: " << interfaces_iter->get_id() << endl;
        (*interfaces_iter).OutputCellid();
        interfaces_iter++;
    }

    vector<Cell>::iterator cells_iter = cells_.begin();
    vector<Cell>::iterator cells_end = cells_.end();
    while (cells_iter != cells_end)
    {
        cout << "cellid: " << cells_iter->get_id() << endl;
        (*cells_iter).OutputInterfaceid();
        cells_iter++;
    }
}