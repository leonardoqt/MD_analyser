#include <cstring>
#include <iomanip>
#include "cell.h"

using namespace std;

void cell :: init(ifstream& in)
{
	string label_tri="xy xz yz";
	string tmp;
	in.clear(); in.seekg(ios::beg);
	// get number of atoms
	getline(in,tmp);
	getline(in,tmp);
	getline(in,tmp);
	in>>num_a;
	getline(in,tmp);
	num_a /= 5;
	num_b = num_a;
	num_c = 3*num_a;
	// check if triclinic or not
	shift.resize(3);
	param.resize(3);
	param[0].resize(3);
	param[1].resize(3);
	param[2].resize(3);
	getline(in,tmp);
	if(tmp.find(label_tri) == string::npos)
	{
		// not triclinic
		in>>shift[0]>>param[0][0];
		getline(in,tmp);
		in>>shift[1]>>param[1][1];
		getline(in,tmp);
		in>>shift[2]>>param[2][2];
		getline(in,tmp);
		param[0][1] = param[0][2] = 0;
		param[1][0] = param[1][2] = 0;
		param[2][0] = param[2][1] = 0;
		param[0][0] -= shift[0];
		param[1][1] -= shift[1];
		param[2][2] -= shift[2];
	}
	else
	{
		// triclinic
		in>>shift[0]>>param[0][0]>>param[1][0];
		getline(in,tmp);
		in>>shift[1]>>param[1][1]>>param[2][0];
		getline(in,tmp);
		in>>shift[2]>>param[2][2]>>param[2][1];
		getline(in,tmp);
		param[0][1] = param[0][2] = 0;
		param[1][2] = 0;
		param[0][0] -= shift[0];
		param[1][1] -= shift[1];
		param[2][2] -= shift[2];
	}
	// allocate all vectors
	A.resize(num_a);
	B.resize(num_b);
	C.resize(num_c);
	oct.resize(num_b);
}

void cell :: read(ifstream& in)
{
	string label_atom = "ITEM: ATOMS";
	string tmp;
	getline(in,tmp);
	while(tmp.find(label_atom) == string::npos)
		getline(in,tmp);
	for(auto& m1 : A)
	{
		in>>m1;
		getline(in,tmp);
	}
	for(auto& m1 : B)
	{
		in>>m1;
		getline(in,tmp);
	}
	for(auto& m1 : C)
	{
		in>>m1;
		getline(in,tmp);
	}
}


void cell :: print()
{
	cout<<"PRIMVEC"<<endl;
	cout<<fixed<<setprecision(6)<<setw(11)<<param[0][0]<<setw(11)<<param[0][1]<<setw(11)<<param[0][2]<<endl;
	cout<<fixed<<setprecision(6)<<setw(11)<<param[1][0]<<setw(11)<<param[1][1]<<setw(11)<<param[1][2]<<endl;
	cout<<fixed<<setprecision(6)<<setw(11)<<param[2][0]<<setw(11)<<param[2][1]<<setw(11)<<param[2][2]<<endl;
	cout<<"PRIMCOORD"<<endl;
	cout<<num_a+num_b+num_c<<endl;
	for(auto& m1 : A) cout<<"Bi "<<m1<<endl;
	for(auto& m1 : B) cout<<"Fe "<<m1<<endl;
	for(auto& m1 : C) cout<<"O  "<<m1<<endl;
}
