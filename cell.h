#ifndef APC523Project_CELL_H
#define APC523Project_CELL_H

#include "model.h"
#include "interface.h"

using namespace std;

class Model;
class Interface;

class Cell
{
public:
	Cell(int *interfaceid, double x, double y, double dx, double dy, int id, Model *mymodel)
	{create_(interfaceid, x, y, dx, dy, id, mymodel);};
	Cell(){create_();};

	void set_cell_interfaces();
		
	int get_id(){return id_;};
	double get_x(){return x_;};
	double get_y(){return y_;};
	double *get_U(){return U_;};
	double get_dt();
	
	void set_U(double*);
	
	void reconstruct(int, char);
	void predict(double, char);
	void update(double, char);
	
	void minbee(double*, double*, double*, double*, double*);
	void superbee(double*, double*, double*, double*, double*);
	
	//for test
	void OutputInterfaceid();

private:
	int interfaceid_[4];
	Interface *cell_interfaces_[4];	// cell_interfaces_[0]: left interface, [1]: right, [2]: bottom, [3]: up
	double x_, y_, dx_, dy_;
	double U_[4];	// U_[0]=rho, U_[1]=rho*u, U_[2]=rho*v, U_[3]=E
	int id_;
	Model *mymodel_;

	void create_(int*, double, double, double, double, int, Model*);
	void create_();
};

#endif
