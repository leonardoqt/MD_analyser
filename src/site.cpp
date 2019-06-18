#include <iomanip>
#include "site.h"

using namespace std;

site :: site()
{
	pos.resize(3);
}

istream& operator>>(istream& in, site&B)
{
	in>>B.pos[0]>>B.pos[1]>>B.pos[2];
	return in;
}
ifstream& operator>>(ifstream& in, site&B)
{
	in>>B.pos[0]>>B.pos[1]>>B.pos[2];
	return in;
}

ostream& operator<<(ostream& out, site&B)
{
	out<<fixed<<setprecision(6)<<setw(11)<<B.pos[0]<<setw(11)<<B.pos[1]<<setw(11)<<B.pos[2];
	return out;
}
ofstream& operator<<(ofstream& out, site&B)
{
	out<<fixed<<setprecision(6)<<setw(11)<<B.pos[0]<<setw(11)<<B.pos[1]<<setw(11)<<B.pos[2];
	return out;
}
