#include <iostream>
#include <fstream>
#include <cstring>
#include "cell.h"

using namespace std;

int main()
{
	ifstream in;
	cell sys1;
	int num_iter;
	string label_iter="ITEM: TIMESTEP";
	string tmp;

	in.open("dump.xyz");
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
	sys1.first_read(in);
	sys1.regist_oct();
	// calculate polarization
	cout<<0<<sys1.ave_p()<<endl;
	for(size_t t1=1; t1<num_iter; t1++)
	{
		sys1.read(in);
		sys1.rebuild_oct();
		cout<<t1<<sys1.ave_p()<<endl;
	}
	return 0;
}
