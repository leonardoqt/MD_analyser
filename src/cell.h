#ifndef __CELL__
#define __CELL__

#include <vector>
#include <iostream>
#include <fstream>
#include "site.h"

using namespace std;

class center
{
public:
	site B;
	vector <site> A;
	vector <site> corr_A;	// = A - site A, usually are 0, but can be non-zero at boundary
	vector <int> ind_A;
	vector <site> C;
	vector <site> corr_C;
	vector <int> ind_C;
};

class cell
{
private:
	vector <double> shift;
	vector < vector <double> > param;
	int num_a;
	int num_b;
	int num_c;
	vector <site> A;
	vector <site> B;
	vector <site> C;
	vector <center> oct;
public:
	void init(ifstream& in);
	void read(ifstream& in);
	void get_oct();		// Bi-Bi < 4.7; Fe-Bi < 4.0; Fe-O < 2.8

	void print();
};
#endif
