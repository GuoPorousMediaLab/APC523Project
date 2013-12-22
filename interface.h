#ifndef APC523Project_INTERFACE_H
#define APC523Project_INTERFACE_H

#include <vector>
#include "cell.h"
#include "model.h"

using namespace std;

class Model;
class Cell;

class Interface
{
public:
	Interface(vector<int> cellid, int id, Model* mymodel){create_(cellid, id, mymodel);};
	Interface(){create_();};

	void set_interface_cells();
	int get_id(){return id_;};

	// for test
	void OutputCellid();

private:
	vector<int> cellid_;
	vector<Cell*> interface_cells_;
	int id_;
	Model* mymodel_;

	void create_(vector<int>, int, Model*);
	void create_();
};


#endif