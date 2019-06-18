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
	const double q_a = 3;
	const double q_b = 3;
	const double q_c = -2;
public:
	void init(ifstream& in);
	void first_read(ifstream& in);
	void read(ifstream& in);
	void regist_oct();		// Bi-Bi < 4.7; Fe-Bi < 4.0; Fe-O < 2.8
	void rebuild_oct();
	site ave_p();

	void print();
};
#endif
