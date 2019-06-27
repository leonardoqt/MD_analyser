#include <iostream>
#include <fstream>
#include <cstring>
#include "cell.h"

using namespace std;

int main()
{
	ifstream in, dat;
	ofstream out_polar, out_power;
	cell sys1;
	int num_iter;
	string label_iter="ITEM: TIMESTEP";
	string tmp;
	site ave_p, std_err_p;

	in.open("dump.xyz");
	dat.open("BFO.data");
	// find number of iterations
	num_iter = 0;
	while(!in.eof())
	{
		getline(in,tmp);
		if(tmp.find(label_iter) != string::npos)
			num_iter++;
	}
	in.clear(); in.seekg(ios::beg);
	// initialize cell
	sys1.init(in);
	sys1.first_read(dat);
	if (!sys1.recover_oct()) sys1.regist_oct();
	// open output file
	out_polar.open("polarization.dat");
#ifdef __SPECTRA__
	sys1.init_spectra(num_iter);
	out_power.open("power_spectra.dat");
#endif
	//===============start calculation=================
	for(size_t t1=0; t1<num_iter; t1++)
	{
		sys1.read(in);
		sys1.rebuild_oct();
		sys1.ave_p(ave_p,std_err_p);
		out_polar<<t1<<ave_p<<std_err_p<<endl;
#ifdef __SPECTRA__
		sys1.save_traj(t1);
#endif
	}
#ifdef __SPECTRA__
	sys1.get_spectra();
	for(size_t t1=0; t1<num_iter/2; t1++)
		out_power<<sys1.freq[t1]<<'\t'<<sys1.spectra[t1]<<endl;
#endif
	return 0;
}
