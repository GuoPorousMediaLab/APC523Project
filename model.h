#ifndef APC523Project_MODEL_H
#define APC523Project_MODEL_H

#include <vector>
#include "cell.h"
#include "interface.h"

using namespace std;

class Cell;
class Interface;

class Model
{
public:
	Model(int Nx, int Ny){create_(Nx, Ny);};
	Model(){create_();};

	vector<Cell> get_cells(){return cells_;};
	vector<Interface> get_interfaces(){return interfaces_;};

	// for test
	void Outputid();

private:
	vector<Cell> cells_;
	vector<Interface> interfaces_;

	int Nx_, Ny_;
	void create_(int, int);
	void create_();

};

#endif