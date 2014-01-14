#include <vector>
#include <cmath>
#include "model.h"

double Model::Timestep(double CPL)
// calculate the step size for next update
{
	double dt = 999.0;
	
    for (vector<int>::size_type i = 0; i < cells_.size(); i++)
    {
		dt = fmin(dt, CPL * cells_[i].get_dt());
	}
	
	return dt;
}

void Model::Reconstructx(int slope_limiter)
// reconstruct each cell in x direction
{
    for (vector<int>::size_type i = 0; i < cells_.size(); i++)
    {
		cells_[i].reconstruct(slope_limiter, 'x');
	}	
}

void Model::Reconstructy(int slope_limiter)
// reconstruct each cell in y direction
{
    for (vector<int>::size_type i = 0; i < cells_.size(); i++)
    {
		cells_[i].reconstruct(slope_limiter, 'y');
	}	
}

void Model::Predictx(double dt)
{
    for (vector<int>::size_type i = 0; i < cells_.size(); i++)
    {
    	cells_[i].predict(dt, 'x');
    }
}

void Model::Predicty(double dt)
{
    for (vector<int>::size_type i = 0; i < cells_.size(); i++)
    {
    	cells_[i].predict(dt, 'y');
    }
}

void Model::Riemannx(int riemann_solver)
// solve riemann problem of each vertical interface
{
	for (int i = 0; i < Nvertical_; i++)
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
// solve riemann problem of each horizontal interface
{
	for (vector<int>::size_type i = Nvertical_; i < interfaces_.size(); i++)
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
    for (vector<int>::size_type i = 0; i < cells_.size(); i++)
	{
		cells_[i].update(dt, 'x');
	}
	
}

void Model::Updatey(double dt)
// update the value of each cell
{
    for (vector<int>::size_type i = 0; i < cells_.size(); i++)
	{
		cells_[i].update(dt, 'y');
	}
	
}
