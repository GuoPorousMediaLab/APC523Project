#include <vector>
#include <cmath>
#include <iostream>
#include <sstream>
#include <fstream>
#include "model.h"

using namespace std;

void Model::create_(const char *filename)
{
	int cellid[2], interfaceid[4];
    double x, y, dx, dy;
    int i;
    string line;
    istringstream data;
    ifstream fin;
    
    fin.open(filename);
    
    getline(fin, line);
    fin >> Nvertical_ >> Nx_ >> Ny_;
    
	getline(fin, line);
    getline(fin, line);

    getline(fin, line);
    while (line[0] != '#')
	{
    	data.str(line);
    	data >> interfaceid[0] >> interfaceid[1] >> interfaceid[2] >> interfaceid[3] >> x >> y >> dx >> dy >> cellid[0];
		cells_.push_back(Cell(interfaceid, x, y, dx, dy, cellid[0], this));
    	//cout << interfaceid[0] << interfaceid[1] << interfaceid[2] << interfaceid[3] << x << y << dx << dy << cellid[0] << endl;
    	getline(fin, line);
    }

    getline(fin, line);
    while (!fin.eof())
	{
		data.str(line);
    	data >> cellid[0] >> cellid[1] >> interfaceid[0];
    	interfaces_.push_back(Interface(cellid, interfaceid[0], this));
    	//cout << cellid[0] << cellid[1] << interfaceid[0] << endl;
    	getline(fin, line);
    }

    // setup cells correspond to interface and interfaces correspond to cells
    for (i = 0; i < cells_.size(); i++)
    {
        cells_[i].set_cell_interfaces();
    }

    for (i = 0; i < interfaces_.size(); i++)
    {
        interfaces_[i].set_interface_cells();
    }

}

void Model::Initialize()
{
// set up a Riemann problem in x direction for testing
	double rho = 1.4, u = 3.0, v = 0.0, p = 1.0;
	double U[4];
	int i;
	
	U[0] = rho;
	U[1] = rho * u;
	U[2] = rho * v;
	U[3] = p / (GAMMA - 1.0) + 0.5 * rho * (u * u + v * v);
	
    for (i = 0; i< cells_.size(); i++)
    {
		cells_[i].initialize(U);
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
//    for (i = 0; i< cells_.size(); i++)
//    {
//    	if (cells_[i].get_y() < y0)
//    	{
//			cells_[i].initialize(Ul);
//    	}
//    	else
//    	{
//			cells_[i].initialize(Ur);
//    	}
//    }

// Set up the flux at fixed-flux boundaries
	for (i = 0; i < Nvertical_; i++)
	{
		interfaces_[i].initialize('x');
	}
	for (i = Nvertical_; i < interfaces_.size(); i++)
	{
		interfaces_[i].initialize('y');
	}
}

double Model::Timestep(double CPL)
// calculate the step size for next update
{
	double dt = 999.0;
	int i;
	
    for (i = 0; i < cells_.size(); i++)
    {
		dt = fmin(dt, CPL * cells_[i].get_dt());
	}
	
	return dt;
}

void Model::Reconstructx(int slope_limiter)
{
	int i;
	
    for (i = 0; i < cells_.size(); i++)
    {
		cells_[i].reconstruct(slope_limiter, 'x');
	}	
}

void Model::Reconstructy(int slope_limiter)
{
	int i;
	
    for (i = 0; i < cells_.size(); i++)
    {
		cells_[i].reconstruct(slope_limiter, 'y');
	}	
}

void Model::Predictx(double dt)
{
	int i;
	
    for (i = 0; i < cells_.size(); i++)
    {
    	cells_[i].predict(dt, 'x');
    }
}

void Model::Predicty(double dt)
{
	int i;
	
    for (i = 0; i < cells_.size(); i++)
    {
    	cells_[i].predict(dt, 'y');
    }
}

void Model::Riemannx(int riemann_solver)
{
	int i;
	
	for (i = 0; i < Nvertical_; i++)
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
	int i;
	
	for (i = Nvertical_; i < interfaces_.size(); i++)
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
	int i;
	
    for (i = 0; i < cells_.size(); i++)
	{
		cells_[i].update(dt, 'x');
	}	
}

void Model::Updatey(double dt)
// update the value of each cell
{
	int i;
	
    for (i = 0; i < cells_.size(); i++)
	{
		cells_[i].update(dt, 'y');
	}	
}

void Model::Outputid()
{
	int i;
	
    for (i = 0; i < interfaces_.size(); i++)
    {
        cout << "interfaceid: " << interfaces_[i].get_id() << endl;
        interfaces_[i].OutputCellid();
    }
    
	cout << endl;
	
    for (i = 0; i < cells_.size(); i++)
    {
        cout << "cellid: " << cells_[i].get_id() << endl;
        cells_[i].OutputInterfaceid();
    }
}

void Model::Outputvalue(const char *filename)
// save the results in a file
{
	double x, y, p, u, v, a, rho;
	double *U;
	int i, j, N;
	ofstream fout;
	
	fout.open(filename);
	fout << "Nx\tNy" << endl;
	fout << Nx_ << '\t' << Ny_ << endl;
	
	fout << "x\ty\tp\tu\tv\ta\trho" << endl;
	N = 0;
	i = 0;
    for (i = 0; i < cells_.size(); i++)
    {
		if (i < cells_.size() - 1 && cells_[i + 1].get_y() != cells_[i].get_y())
		{
			for (j = N; j < Nx_; j++)
			{
				x = cells_[i - 1].get_x() + (j - N + 1) * cells_[i - 1].get_dx();
				y = cells_[i - 1].get_y();
				rho = 0.0;
				p = 0.0;
				u = 0.0;
				v = 0.0;
				a = 0.0;
				fout << x << '\t' << y << '\t' << p << '\t' << u << '\t' << v << '\t' << a << '\t' << rho << endl;
			}
			N = 0;
		}
		else
		{
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
		}
	}
	
	fout.close();
	
	double x1, x2, y1, y2;
	fout.open("interfaces.txt");
	for (i = 0; i < Nvertical_; i++)
	{
		if (interfaces_[i].get_cellid(0) < 0)
		{
			x1 = interfaces_[i].get_cell(1)->get_x() - interfaces_[i].get_cell(1)->get_dx() / 2.0;
			x2 = x1;
			y1 = interfaces_[i].get_cell(1)->get_y() - interfaces_[i].get_cell(1)->get_dy() / 2.0;
			y2 = interfaces_[i].get_cell(1)->get_y() + interfaces_[i].get_cell(1)->get_dy() / 2.0;
			fout << x1 << '\t' << x2 << '\t' << y1 << '\t' << y2 << '\t' << interfaces_[i].get_cellid(0) << endl;
		}
		else if (interfaces_[i].get_cellid(1) < 0)
		{
			x1 = interfaces_[i].get_cell(0)->get_x() + interfaces_[i].get_cell(0)->get_dx() / 2.0;
			x2 = x1;
			y1 = interfaces_[i].get_cell(0)->get_y() - interfaces_[i].get_cell(0)->get_dy() / 2.0;
			y2 = interfaces_[i].get_cell(0)->get_y() + interfaces_[i].get_cell(0)->get_dy() / 2.0;
			fout << x1 << '\t' << x2 << '\t' << y1 << '\t' << y2 << '\t' << interfaces_[i].get_cellid(1) << endl;
		}
		else
		{
			x1 = interfaces_[i].get_cell(0)->get_x() + interfaces_[i].get_cell(0)->get_dx() / 2.0;
			x2 = interfaces_[i].get_cell(1)->get_x() - interfaces_[i].get_cell(1)->get_dx() / 2.0;
			y1 = interfaces_[i].get_cell(0)->get_y() - interfaces_[i].get_cell(0)->get_dy() / 2.0;
			y2 = interfaces_[i].get_cell(1)->get_y() + interfaces_[i].get_cell(1)->get_dy() / 2.0;
			fout << x1 << '\t' << x2 << '\t' << y1 << '\t' << y2 << '\t' << interfaces_[i].get_cellid(0) << endl;
		}
	}
	for (i = Nvertical_; i < interfaces_.size(); i++)
	{
		if (interfaces_[i].get_cellid(0) < 0)
		{
			x1 = interfaces_[i].get_cell(1)->get_x() - interfaces_[i].get_cell(1)->get_dx() / 2.0;
			x2 = interfaces_[i].get_cell(1)->get_x() + interfaces_[i].get_cell(1)->get_dx() / 2.0;
			y1 = interfaces_[i].get_cell(1)->get_y() - interfaces_[i].get_cell(1)->get_dy() / 2.0;
			y2 = y1;
			fout << x1 << '\t' << x2 << '\t' << y1 << '\t' << y2 << '\t' << interfaces_[i].get_cellid(0) << endl;
		}
		else if (interfaces_[i].get_cellid(1) < 0)
		{
			x1 = interfaces_[i].get_cell(0)->get_x() - interfaces_[i].get_cell(0)->get_dx() / 2.0;
			x2 = interfaces_[i].get_cell(0)->get_x() + interfaces_[i].get_cell(0)->get_dx() / 2.0;
			y1 = interfaces_[i].get_cell(0)->get_y() + interfaces_[i].get_cell(0)->get_dy() / 2.0;
			y2 = y1;
			fout << x1 << '\t' << x2 << '\t' << y1 << '\t' << y2 << '\t' << interfaces_[i].get_cellid(1) << endl;
		}
		else
		{
			x1 = interfaces_[i].get_cell(0)->get_x() - interfaces_[i].get_cell(0)->get_dx() / 2.0;
			x2 = interfaces_[i].get_cell(1)->get_x() + interfaces_[i].get_cell(1)->get_dx() / 2.0;
			y1 = interfaces_[i].get_cell(0)->get_y() + interfaces_[i].get_cell(0)->get_dy()/ 2.0;
			y2 = interfaces_[i].get_cell(1)->get_y() - interfaces_[i].get_cell(1)->get_dy()/ 2.0;;
			fout << x1 << '\t' << x2 << '\t' << y1 << '\t' << y2 << '\t' << interfaces_[i].get_cellid(0) << endl;
		}
	}
	
	
}
