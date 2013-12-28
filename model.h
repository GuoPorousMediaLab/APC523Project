#ifndef APC523Project_MODEL_H
#define APC523Project_MODEL_H

#include <vector>
#include "cell.h"
#include "interface.h"

#define GAMMA 1.4

using namespace std;

class Cell;
class Interface;

class Model
{
public:
	Model(int Nx, int Ny){create_(Nx, Ny);};
	Model(){create_();};

	Cell *get_cell(int i){return &cells_[i];};
	Interface *get_interface(int i){return &interfaces_[i];};
	
	void Initialize();
	void Reconstructx(int);
	void Riemannx(int);
	void Updatex(double);
//	void Reorder();
//	void Reconstructy(int);
//	void Riemanny(int);
//	void Updatey(double);
	double Timestep(double);
	
	// for test
	void Outputid();
	void Outputvalue(const char*);
	

private:
	int Nx_, Ny_;
	vector<Cell> cells_;
	vector<Interface> interfaces_;
	void create_(int, int);
	void create_();
};

#endif
