#include <iomanip>
#include "site.h"

using namespace std;

site :: site()
{
	pos.resize(3);
}

void site :: clear()
{
	pos[0] = pos[1] = pos[2] = 0;
}

site site :: operator+(const site&B)
{
	site res;
	res.pos[0] = pos[0] + B.pos[0];
	res.pos[1] = pos[1] + B.pos[1];
	res.pos[2] = pos[2] + B.pos[2];
	return res;
}

site site :: operator-(const site&B)
{
	site res;
	res.pos[0] = pos[0] - B.pos[0];
	res.pos[1] = pos[1] - B.pos[1];
	res.pos[2] = pos[2] - B.pos[2];
	return res;
}

site site :: operator*(const double&B)
{
	site res;
	res.pos[0] = pos[0]*B;
	res.pos[1] = pos[1]*B;
	res.pos[2] = pos[2]*B;
	return res;
}

site site :: operator/(const double&B)
{
	site res;
	res.pos[0] = pos[0]/B;
	res.pos[1] = pos[1]/B;
	res.pos[2] = pos[2]/B;
	return res;
}

istream& operator>>(istream& in, site& B)
{
	in>>B.pos[0]>>B.pos[1]>>B.pos[2];
	return in;
}
ifstream& operator>>(ifstream& in, site& B)
{
	in>>B.pos[0]>>B.pos[1]>>B.pos[2];
	return in;
}

ostream& operator<<(ostream& out, site B)
{
	out<<fixed<<setprecision(6)<<setw(11)<<B.pos[0]<<setw(11)<<B.pos[1]<<setw(11)<<B.pos[2];
	return out;
}
ofstream& operator<<(ofstream& out, site B)
{
	out<<fixed<<setprecision(6)<<setw(11)<<B.pos[0]<<setw(11)<<B.pos[1]<<setw(11)<<B.pos[2];
	return out;
}
