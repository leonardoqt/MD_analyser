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
	//	For BFO, Bi=4.68 Fe=3.66 O=-2.78
	const double q_a = 4.68;
	const double q_b = 3.66;
	const double q_c = -2.78;
#ifdef __SPECTRA__
	double dw;
	double **traj;
#endif

public:
	const double dt = 0.04;
	int tot_step;
	vector <site> polarization;
	vector <site> polarization_stderr;
	void init(ifstream& in);
	void first_read(ifstream& in);
	void read(ifstream& in);
	void regist_oct();		// Bi-Bi < 4.7; Fe-Bi < 4.0; Fe-O < 2.8
	int  recover_oct();
	void rebuild_oct();
	void ave_p(int iter);
#ifdef __SPECTRA__
	vector <double> freq;
	vector <site> p_w;	//	fft of polarization
	vector <double> spectra;
	void get_p_w();
	void init_spectra();
	void save_traj(int iter);
	void get_spectra();
#endif

	void print();
};
#endif
