#ifndef APC523Project_CELL_H
#define APC523Project_CELL_H

#include <vector>
#include "model.h"
#include "interface.h"

using namespace std;

class Model;
class Interface;

class Cell
{
public:
	Cell(vector<int> interface_id, double x, double y, double dx, double dy, int id, Model* mymodel)
	{create_(interface_id, x, y, dx, dy, id, mymodel);};
	Cell(){create_();};

	void set_cell_interfaces();
	int get_id(){return id_;};

	void initialize();

	//for test
	void OutputInterfaceid();

private:
	vector<int> interfaceid_;
	vector<Interface*> cell_interfaces_;
	double x_, y_, dx_, dy_;
	int id_;
	Model* mymodel_;

	vector<double> U_;

	void create_(vector<int>, double, double, double, double, int, Model*);
	void create_();

};

#endif