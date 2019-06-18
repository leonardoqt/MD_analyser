#include <iostream>
#include <fstream>
#include <cstring>
#include "cell.h"

using namespace std;

int main()
{
	ifstream in, dat;
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
	sys1.regist_oct();
	// calculate polarization
	for(size_t t1=0; t1<num_iter; t1++)
	{
		sys1.read(in);
		sys1.rebuild_oct();
		sys1.ave_p(ave_p,std_err_p);
		cout<<t1<<ave_p<<std_err_p<<endl;
	}
	return 0;
}
