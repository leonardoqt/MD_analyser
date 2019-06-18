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
	sys1.rebuild_oct();
//	sys1.print();
	return 0;
}
