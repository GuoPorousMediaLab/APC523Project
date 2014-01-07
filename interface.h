#ifndef APC523Project_INTERFACE_H
#define APC523Project_INTERFACE_H

#include "cell.h"
#include "model.h"

using namespace std;

class Model;
class Cell;

class Interface
{
public:
	Interface(int *cellid, int id, Model *mymodel){create_(cellid, id, mymodel);};
	Interface(){create_();};

	void set_interface_cells();
	void initialize();
	
	int get_id(){return id_;};
	Cell *get_cell(int i){return interface_cells_[i];};
	int get_cellid(int i){return cellid_[i];}
	double *get_F(){return F_;};
	double *get_U1(){return U1_;};
	double *get_U2(){return U2_;};
	
	void set_U1(double *);
	void set_U2(double *);
	void roe(char);
	void hlle(char);
	void hllc(char);
	
	// for test
	void OutputCellid();

private:
	int cellid_[2];
	Cell *interface_cells_[2];	// interface_cells_[0]: bottom/left cell, [1] : top/right cell
	int id_;
	Model *mymodel_;

	// U1_ correspondes to cell1 of the interface, U2_ corresponds to cell2
	double U1_[4];
	double U2_[4];
	// direction of flux F is from cell1 to cell2
	double F_[4];

	void create_(int*, int, Model*);
	void create_();
};


#endif
