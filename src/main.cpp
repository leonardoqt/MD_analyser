#include <iostream>
#include <iomanip>
#include <fstream>
#include <cstring>
#include "cell.h"

using namespace std;

int main(int argc, char* argv[])
{
	ifstream in, dat;
	ofstream out_polar, out_polar_fft, out_power;
	cell sys1;
	string tmp;
	site ave_p, std_err_p;

	//in.open("dump.xyz");
	//dat.open("BFO.data");
	in.open(argv[1]);
	dat.open(argv[2]);
	// initialize cell
	sys1.init(in);
	sys1.first_read(dat);
	if (!sys1.recover_oct()) sys1.regist_oct();
	// open output file
	out_polar.open("polarization.dat");
#ifdef __SPECTRA__
	sys1.init_spectra();
	out_polar_fft.open("fft_polarization.dat");
	out_power.open("power_spectra.dat");
#endif
	//===============start calculation=================
	//cout<<"Start calculating polarization..."<<endl<<"[0/"<<sys1.tot_step<<']';
	cout<<"Start calculating polarization..."<<endl<<"[0/"<<sys1.tot_step<<']'<<endl;
	for(size_t t1=0; t1<sys1.tot_step; t1++)
	{
		sys1.read(in);
		sys1.rebuild_oct();
		sys1.ave_p(t1);
#ifdef __SPECTRA__
		sys1.save_traj(t1);
#endif
		//cout<<'\r'<<'['<<t1+1<<'/'<<sys1.tot_step<<']';
		cout<<"\e[A"<<'['<<t1+1<<'/'<<sys1.tot_step<<']'<<endl;
	}
	for(size_t t1=0; t1<sys1.tot_step; t1++)
		out_polar<<fixed<<setprecision(3)<<setw(7)<<t1*sys1.dt<<sys1.polarization[t1]<<sys1.polarization_stderr[t1]<<endl;
#ifdef __SPECTRA__
	sys1.get_p_w();
	sys1.get_spectra();
	for(size_t t1=0; t1<sys1.tot_step/2; t1++)
		out_polar_fft<<fixed<<setprecision(3)<<setw(7)<<sys1.freq[t1]<<sys1.p_w[t1]<<endl;
	for(size_t t1=0; t1<sys1.tot_step/2; t1++)
		out_power<<fixed<<setprecision(3)<<setw(7)<<sys1.freq[t1]<<setprecision(6)<<setw(13)<<sys1.spectra[t1]<<endl;
#endif
	cout<<"Job done"<<endl;
	return 0;
}
