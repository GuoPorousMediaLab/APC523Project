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
	double centerx = 0.0;
    double centery = 0.0;
    double dx = 1.0 / Nx;
    double dy = 1.0 / Ny;
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
                // -1 -- fixed boundary
                // -2 -- transimissive boundary
                // -3 -- reflective boundary
            }
            if (j < Nx_)
			{
                cellid[1] = i * Nx_ + j;	// right cell
            }
            else
            {
            	cellid[1] = -2;	// interface is right boundary
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
            	cellid[0] = -2;	// interface is bottom boundary
            }
            if (i < Ny_)
			{
                cellid[1] = i * Nx_ + j;	// top cell
            }
            else
            {
            	cellid[1] = -2;	// interface is top boundary
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
{
// set up a Riemann problem in x direction for testing
	double rho1 = 1.4, u1 = 2.0, v1 = 0.0, p1 = 1.0, rho2 = 1.4, u2 = 0.0, v2 = 0.0, p2 = 1.0, x0 = 0.0;
	double Ul[4], Ur[4];
	int i;
	
	Ul[0] = rho1;
	Ul[1] = rho1 * u1;
	Ul[2] = rho1 * v1;
	Ul[3] = p1 / (GAMMA - 1.0) + 0.5 * rho1 * (u1 * u1 + v1 * v1);
	Ur[0] = rho2;
	Ur[1] = rho2 * u2;
	Ur[2] = rho2 * v2;
	Ur[3] = p2 / (GAMMA - 1.0) + 0.5 * rho2 * (u2 * u2 + v2 * v2);
	
	vector<Cell>::iterator cells_iter;
    for (cells_iter = cells_.begin(); cells_iter != cells_.end(); cells_iter++)
    {
    	if (cells_iter->get_x() < x0)
    	{
			cells_iter->initialize(Ul);
    	}
    	else
    	{
			cells_iter->initialize(Ur);
    	}
    }

// set up a Riemann problem in y direction for testing
//	double rho1 = 1.0, u1 = 0.1, v1 = 0.75, p1 = 1.0, rho2 = 0.125, u2 = 0.5, v2 = 0.0, p2 = 0.1, y0 = -0.2;
//	double Ul[4], Ur[4];
//	int i;
//	
//	Ul[0] = rho1;
//	Ul[1] = rho1 * u1;
//	Ul[2] = rho1 * v1;
//	Ul[3] = p1 / (GAMMA - 1.0) + 0.5 * rho1 * (u1 * u1 + v1 * v1);
//	Ur[0] = rho2;
//	Ur[1] = rho2 * u2;
//	Ur[2] = rho2 * v2;
//	Ur[3] = p2 / (GAMMA - 1.0) + 0.5 * rho2 * (u2 * u2 + v2 * v2);
//	
//	vector<Cell>::iterator cells_iter;
//    for (cells_iter = cells_.begin(); cells_iter != cells_.end(); cells_iter++)
//    {
//    	if (cells_iter->get_y() < y0)
//    	{
//			cells_iter->initialize(Ul);
//    	}
//    	else
//    	{
//			cells_iter->initialize(Ur);
//    	}
//    }

// set up a "tilted" Riemann problem for 2d testing
//	double rho1 = 1.0, u1 = 0.75 / sqrt(2), v1 = -0.75 / sqrt(2), p1 = 1.0, rho2 = 0.125, u2 = 0.0, v2 = 0.0, p2 = 0.1;
//	double Ul[4], Ur[4];
//	int i;
//	
//	Ul[0] = rho1;
//	Ul[1] = rho1 * u1;
//	Ul[2] = rho1 * v1;
//	Ul[3] = p1 / (GAMMA - 1.0) + 0.5 * rho1 * (u1 * u1 + v1 * v1);
//	Ur[0] = rho2;
//	Ur[1] = rho2 * u2;
//	Ur[2] = rho2 * v2;
//	Ur[3] = p2 / (GAMMA - 1.0) + 0.5 * rho2 * (u2 * u2 + v2 * v2);
//	
//	vector<Cell>::iterator cells_iter;
//    for (cells_iter = cells_.begin(); cells_iter != cells_.end(); cells_iter++)
//    {
//    	if (cells_iter->get_y() - cells_iter->get_x() - 0.2 * sqrt(2) > 0.0)
//    	{
//			cells_iter->initialize(Ul);
//    	}
//    	else
//    	{
//			cells_iter->initialize(Ur);
//    	}
//    }

// Set up the flux at fixed-flux boundaries
	for (int i = 0; i < (Nx_ + 1) * Ny_; i++)
	{
		interfaces_[i].initialize('x');
	}
	for (int i = (Nx_ + 1) * Ny_; i < 2 * Nx_ * Ny_ + Nx_ + Ny_; i++)
	{
		interfaces_[i].initialize('y');
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

void Model::Reconstructx(int slope_limiter)
{
	vector<Cell>::iterator cells_iter;
	
    for (cells_iter = cells_.begin(); cells_iter != cells_.end(); cells_iter++)
    {
		cells_iter->reconstruct(slope_limiter, 'x');
	}	
}

void Model::Reconstructy(int slope_limiter)
{
	vector<Cell>::iterator cells_iter;
	
    for (cells_iter = cells_.begin(); cells_iter != cells_.end(); cells_iter++)
    {
		cells_iter->reconstruct(slope_limiter, 'y');
	}	
}

void Model::Predictx(double dt)
{
	vector<Cell>::iterator cells_iter;
	
    for (cells_iter = cells_.begin(); cells_iter != cells_.end(); cells_iter++)
    {
    	cells_iter->predict(dt, 'x');
    }
}

void Model::Predicty(double dt)
{
	vector<Cell>::iterator cells_iter;
	
    for (cells_iter = cells_.begin(); cells_iter != cells_.end(); cells_iter++)
    {
    	cells_iter->predict(dt, 'y');
    }
}

void Model::Riemannx(int riemann_solver)
{
	for (int i = 0; i < (Nx_ + 1) * Ny_; i++)
	{
		switch (riemann_solver)
		{
			case 1:
				interfaces_[i].roe('x');
				break;
			case 2:
				interfaces_[i].hlle('x');
				break;
			case 3:
				interfaces_[i].hllc('x');
				break;
			default:
				interfaces_[i].roe('x');
		}
	}	
}

void Model::Riemanny(int riemann_solver)
{
	for (int i = (Nx_ + 1) * Ny_; i < 2 * Nx_ * Ny_ + Nx_ + Ny_; i++)
	{
		switch (riemann_solver)
		{
			case 1:
				interfaces_[i].roe('y');
				break;
			case 2:
				interfaces_[i].hlle('y');
				break;
			case 3:
				interfaces_[i].hllc('y');
				break;
			default:
				interfaces_[i].roe('y');
		}
	}	
}

void Model::Updatex(double dt)
// update the value of each cell
{
	vector<Cell>::iterator cells_iter;
	
    for (cells_iter = cells_.begin(); cells_iter != cells_.end(); cells_iter++)
	{
		cells_iter->update(dt, 'x');
	}	
}

void Model::Updatey(double dt)
// update the value of each cell
{
	vector<Cell>::iterator cells_iter;
	
    for (cells_iter = cells_.begin(); cells_iter != cells_.end(); cells_iter++)
	{
		cells_iter->update(dt, 'y');
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
