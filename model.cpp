#include <vector>
#include <cmath>
#include <iostream>
#include <fstream>
#include "model.h"

using namespace std;

void Model::create_(int Nx, int Ny)
{
	Nx_ = Nx;
	Ny_ = Ny;

	// initialize cells
	int interfaceid[4]; // interfaceid for each cell
	double centerx = 0.5;
    double centery = 0.0;
    double dx = 0.01;
    double dy = 0.01;
    double startx = centerx - (Nx_ - 1) / 2.0 * dx;
    double starty = centery - (Ny_ - 1) / 2.0 * dy;
    int i, j;
    for (i = 0; i < Ny_; i++)
    {    
    	for (j = 0; j < Nx_; j++)
        {
        	interfaceid[0] = i * (Nx_ + 1) + j;	// left interface
        	interfaceid[1] = i * (Nx_ + 1) + j + 1;	// right interface
        	interfaceid[2] =  Ny_ * (Nx_ + 1) + i * Nx_ + j;	// bottom interface
        	interfaceid[3] =  Ny_ * (Nx_ + 1) + (i + 1) * Nx_ + j;	// up interface
            cells_.push_back(Cell(interfaceid, startx + j * dx, starty + i * dy, dx, dy, i * Nx_ + j, this));
        }
    }

	// initialize intefaces
	int cellid[2]; // cellid for each interface
	
	// vertical interfaces
	for (int i = 0; i < Ny_; i++)
	{
		for (int j = 0; j < Nx_ + 1; j++)
		{
			if (j > 0)
			{
                cellid[0] = i * Nx_ + j - 1;	// left cell
            }
            else
            {
            	cellid[0] = -1;	// interface is left boundary
            }
            if (j < Nx_)
			{
                cellid[1] = i * Nx_ + j;	// right cell
            }
            else
            {
            	cellid[1] = -1;	// interface is right boundary
            }
			interfaces_.push_back(Interface(cellid, i * (Nx_ + 1) + j, this));
		}
	}

	// horizontal interfaces
    for (int i = 0; i < Ny_ + 1; i++)
    {
        for (int j = 0; j < Nx_; j++)
        {
            if (i > 0)
			{
                cellid[0] = (i - 1) * Nx_ + j;	// bottom cell
            }
            else
            {
            	cellid[0] = -1;	// interface is bottom boundary
            }
            if (i < Ny_)
			{
                cellid[1] = i * Nx_ + j;	// top cell
            }
            else
            {
            	cellid[1] = -1;	// interface is top boundary
            }
            interfaces_.push_back(Interface(cellid, (Nx_ + 1) * Ny_ + i * Nx_ + j, this));
        }
    }

    // setup cells correspond to interface and interfaces correspond to cells
    vector<Interface>::iterator interfaces_iter;
    for (interfaces_iter = interfaces_.begin(); interfaces_iter != interfaces_.end(); interfaces_iter++)
    {
        interfaces_iter->set_interface_cells();
    }

    vector<Cell>::iterator cells_iter;
    for (cells_iter = cells_.begin(); cells_iter != cells_.end(); cells_iter++)
    {
        cells_iter->set_cell_interfaces();
    }

}

void Model::create_()
{
	Nx_ = 0;
	Ny_ = 0;
}

void Model::Initialize()
// set up a Riemann problem for testing
// U = Ul when x < x0 and U = Ur when x > x0
{
	double rho1 = 1.0, u1 = 0.75, p1 = 1.0, rho2 = 0.125, u2 = 0.0, p2 = 0.1, x0 = 0.3;
	double Ul[4], Ur[4];
	int i;
	
	Ul[0] = rho1;
	Ul[1] = rho1 * u1;
	Ul[2] = 0.0;
	Ul[3] = p1 / (GAMMA - 1.0) + 0.5 * rho1 * u1 * u1;
	Ur[0] = rho2;
	Ur[1] = rho2 * u2;
	Ur[2] = 0.0;
	Ur[3] = p2 / (GAMMA - 1.0) + 0.5 * rho2 * u2 * u2;
	
	vector<Cell>::iterator cells_iter;
    for (cells_iter = cells_.begin(); cells_iter != cells_.end(); cells_iter++)
    {
    	if (cells_iter->get_x() < x0)
    	{
			cells_iter->set_U(Ul);
    	}
    	else
    	{
			cells_iter->set_U(Ur);
    	}
    }
}

double Model::Timestep(double CPL)
// calculate the step size for next update
{
	double dt = 999.0;
	vector<Cell>::iterator cells_iter;
	
    for (cells_iter = cells_.begin(); cells_iter != cells_.end(); cells_iter++)
    {
		dt = fmin(dt, CPL * cells_iter->get_dt());
	}
	
	return dt;
}

void Model::Reconstructx(int slopeLimiter)
{
	vector<Cell>::iterator cells_iter;
	
    for (cells_iter = cells_.begin(); cells_iter != cells_.end(); cells_iter++)
    {
//		switch (slopeLimiter)
//		{
//			case 1:
//				cells_iter->minbeex();
//				break;
//			case 2:
//				cells_iter->superbeex();
//				break;
//			default:
//				cells_iter->minbeex();
//		}
		//cells_iter->minbeex();
		cells_iter->superbeex();
	}	
}

void Model::Predictx(double dt)
{
	vector<Cell>::iterator cells_iter;
	
    for (cells_iter = cells_.begin(); cells_iter != cells_.end(); cells_iter++)
    {
    	cells_iter->predictx(dt);
    }
}

void Model::Riemannx(int riemannSolver)
{
	for (int i = 0; i < (Nx_ + 1) * Ny_; i++)
	{
//		switch (riemannSolver)
//		{
//			case 1:
//				interfaces_[i].roe();
//				break;
//			case 2:
//				interfaces_[i].hlle();
//				break;
//			case 3:
//				interfaces_[i].hllc();
//				break;
//			default:
//				interfaces_[i].roe();
//		}
		//interfaces_[i].roe();
		//interfaces_[i].hlle();
		interfaces_[i].hllc();
	}	
}

void Model::Updatex(double dt)
// update the value of each cell
{
	vector<Cell>::iterator cells_iter;
	
    for (cells_iter = cells_.begin(); cells_iter != cells_.end(); cells_iter++)
	{
		cells_iter->updatex(dt);
	}	
}

void Model::Outputid()
{
    vector<Interface>::iterator interfaces_iter;
    for (interfaces_iter = interfaces_.begin(); interfaces_iter != interfaces_.end(); interfaces_iter++)
    {
        cout << "interfaceid: " << interfaces_iter->get_id() << endl;
        interfaces_iter->OutputCellid();
    }
    
	cout << endl;
	
    vector<Cell>::iterator cells_iter;
    for (cells_iter = cells_.begin(); cells_iter != cells_.end(); cells_iter++)
    {
        cout << "cellid: " << cells_iter->get_id() << endl;
        cells_iter->OutputInterfaceid();
    }
}

void Model::Outputvalue(const char *filename)
// save the results in a file
{
	double p, u, v, a, rho;
	double *U;
	int i;
	ofstream fs;
	
	fs.open(filename);
	fs << "x\ty\tp\tu\tv\ta\trho" << endl;
	
    vector<Cell>::iterator cells_iter;
    for (cells_iter = cells_.begin(); cells_iter != cells_.end(); cells_iter++)
    {
		U = cells_iter->get_U();
		
		rho = U[0];
		p = (GAMMA - 1.0) * (U[3] - 0.5 * (U[1] * U[1] + U[2] * U[2])/ U[0]);
		u = U[1] / U[0];
		v = U[2] / U[0];
		a = sqrt(GAMMA * p / U[0]);
		
		fs << cells_iter->get_x() << "\t" << cells_iter->get_y() << "\t" << p << "\t" << u << "\t" << v << "\t" << a << "\t" << rho << endl;
	}
	
	fs.close();
}
