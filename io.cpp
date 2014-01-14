#include <vector>
#include <cmath>
#include <iostream>
#include <sstream>
#include <fstream>
#include <cstdlib>
#include "model.h"

using namespace std;

void Model::create_(const char *filename)
// create the grid
{
	int cellid[2], interfaceid[4];
    double x, y, dx, dy;
    string line;
    ifstream fin;
    
    // read grid information from file
    fin.open(filename);
    
    if (fin.is_open())
    {
	    getline(fin, line);
	    fin >> Nvertical_ >> Nx_ >> Ny_;
	    
		getline(fin, line);
	    getline(fin, line);
		
		// read cell data
	    getline(fin, line);
	    while (line[0] != '#')
		{
			istringstream data(line);
	    	data >> interfaceid[0] >> interfaceid[1] >> interfaceid[2] >> interfaceid[3] >> x >> y >> dx >> dy >> cellid[0];
			cells_.push_back(Cell(interfaceid, x, y, dx, dy, cellid[0], this));
	    	getline(fin, line);
	    }
	
		// read interface data
	    getline(fin, line);
	    while (!fin.eof())
		{
			istringstream data(line);
	    	data >> cellid[0] >> cellid[1] >> interfaceid[0];
	    	interfaces_.push_back(Interface(cellid, interfaceid[0], this));
	    	getline(fin, line);
	    }
	    fin.close();
	}
	else
	{
		cerr << "Can't open " << filename << '.' << endl;
		exit(1);
	}
    // set cell's neighboring interface and interface's neighboring cells
    for (vector<int>::size_type i = 0; i < cells_.size(); i++)
    {
        cells_[i].set_cell_interfaces();
    }

    for (vector<int>::size_type i = 0; i < interfaces_.size(); i++)
    {
        interfaces_[i].set_interface_cells();
    }
}

void Model::Initialize(double mach)
// set a uniform initial condition corresponding to the Mach number
{
	double rho = 1.4, u = mach, v = 0.0, p = 1.0;
	double U[4];
	
	U[0] = rho;
	U[1] = rho * u;
	U[2] = rho * v;
	U[3] = p / (GAMMA - 1.0) + 0.5 * rho * (u * u + v * v);
	
    for (vector<int>::size_type i = 0; i< cells_.size(); i++)
    {
		cells_[i].initialize(U);
    }

// Set up the flux for fixed-flux boundaries
	for (int i = 0; i < Nvertical_; i++)
	{
		interfaces_[i].initialize('x');
	}
	for (vector<int>::size_type i = Nvertical_; i < interfaces_.size(); i++)
	{
		interfaces_[i].initialize('y');
	}
}

void Model::Output(const char *filename)
// save the results in a file
{
	double x, y, p, u, v, a, rho;
	double *U;
	int j, N;
	ofstream fout;
	
	fout.open(filename);
	fout << "# Nx\tNy" << endl;
	fout << Nx_ << '\t' << Ny_ << endl;
	
	fout << "# x\ty\tp\tu\tv\ta\trho" << endl;
	N = 0;
    for (vector<int>::size_type i = 0; i < cells_.size(); i++)
    {
		// read x, y, rho, p, u, v, a from each cell and write them into the file
		x = cells_[i].get_x();
		y = cells_[i].get_y();
		U = cells_[i].get_U();
		rho = U[0];
		p = (GAMMA - 1.0) * (U[3] - 0.5 * (U[1] * U[1] + U[2] * U[2])/ U[0]);
		u = U[1] / U[0];
		v = U[2] / U[0];
		a = sqrt(GAMMA * p / U[0]);
		fout << x << '\t' << y << '\t' << p << '\t' << u << '\t' << v << '\t' << a << '\t' << rho << endl;
		N++;
		
		// add the points inside the cylinder to meet the requirement of MATLAB plot functions
		if (i < cells_.size() - 1 && cells_[i + 1].get_y() != cells_[i].get_y())
		{
			for (j = N; j < Nx_; j++)
			{
				x = cells_[i].get_x() + (j - N + 1) * cells_[i].get_dx();
				y = cells_[i].get_y();
				rho = 0.0;
				p = 0.0;
				u = 0.0;
				v = 0.0;
				a = 0.0;
				fout << x << '\t' << y << '\t' << p << '\t' << u << '\t' << v << '\t' << a << '\t' << rho << endl;
			}
			N = 0;
		}
	}
	fout.close();
}
