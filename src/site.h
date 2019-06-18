#ifndef __SITE__
#define __SITE__

#include <iostream>
#include <fstream>
#include <vector>

using namespace std;

class site
{
public:
	vector <double> pos;

	site();
	friend istream& operator>>(istream&, site&);
	friend ostream& operator<<(ostream&, site&);
	friend ifstream& operator>>(ifstream&, site&);
	friend ofstream& operator<<(ofstream&, site&);
};
#endif
