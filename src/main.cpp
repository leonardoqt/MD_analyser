#include <iostream>
#include <fstream>
#include <cstring>
#include "cell.h"

using namespace std;

int main()
{
	ifstream in;
	cell sys1;

	in.open("dump.xyz");
	sys1.init(in);
	sys1.first_read(in);
	sys1.regist_oct();
	for(size_t t1=0; t1<1000; t1++)
		sys1.read(in);
	sys1.rebuild_oct();
//	sys1.print();
	return 0;
}
