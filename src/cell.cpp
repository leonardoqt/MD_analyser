#include <cstring>
#include <iomanip>
#include <cmath>
#include "cell.h"

using namespace std;

void cell :: init(ifstream& in)
{
	string label_num = "ITEM: NUMBER OF ATOMS";
	string tmp;
	// get number of atoms
	getline(in,tmp);
	while(tmp.find(label_num) == string::npos)
		getline(in,tmp);
	in>>num_a;
	num_a /= 5;
	num_b = num_a;
	num_c = 3*num_a;
	// allocate cell
	shift.resize(3);
	param.resize(3);
	param[0].resize(3);
	param[1].resize(3);
	param[2].resize(3);
	// allocate all vectors
	A.resize(num_a);
	B.resize(num_b);
	C.resize(num_c);
	oct.resize(num_b);
	// rewind ifstream
	in.clear(); in.seekg(ios::beg);
}

void cell :: first_read(ifstream& in)
{
	string label_cell = "ITEM: BOX";
	string label_tri = "xy xz yz";
	string label_atom = "ITEM: ATOMS";
	string tmp;
	// find cell parameter
	getline(in,tmp);
	while(tmp.find(label_cell) == string::npos)
		getline(in,tmp);
	// check if triclinic or not
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
	// get position of atoms
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

void cell :: read(ifstream& in)
{
	string label_cell = "ITEM: BOX";
	string label_tri = "xy xz yz";
	string label_atom = "ITEM: ATOMS";
	string tmp;
	site tmp_site, tmp_corr;
	// find cell parameter
	getline(in,tmp);
	while(tmp.find(label_cell) == string::npos)
		getline(in,tmp);
	// check if triclinic or not
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
	// get position of atoms
	getline(in,tmp);
	while(tmp.find(label_atom) == string::npos)
		getline(in,tmp);
	for(auto& m1 : A)
	{
		in>>tmp_site;
		getline(in,tmp);
		if(tmp_site.pos[0] - m1.pos[0] >= param[0][0]/2)
			tmp_corr.pos[0] = -1;
		else if(tmp_site.pos[0] - m1.pos[0] <= -param[0][0]/2)
			tmp_corr.pos[0] = 1;
		else
			tmp_corr.pos[0] = 0;
		if(tmp_site.pos[1] - m1.pos[1] >= param[1][1]/2)
			tmp_corr.pos[1] = -1;
		else if(tmp_site.pos[1] - m1.pos[1] <= -param[1][1]/2)
			tmp_corr.pos[1] = 1;
		else
			tmp_corr.pos[1] = 0;
		if(tmp_site.pos[2] - m1.pos[2] >= param[2][2]/2)
			tmp_corr.pos[2] = -1;
		else if(tmp_site.pos[2] - m1.pos[2] <= -param[2][2]/2)
			tmp_corr.pos[2] = 1;
		else
			tmp_corr.pos[2] = 0;
		m1.pos[0] = tmp_site.pos[0] + tmp_corr.pos[0]*param[0][0]+tmp_corr.pos[1]*param[1][0]+tmp_corr.pos[2]*param[2][0];
		m1.pos[1] = tmp_site.pos[1] + tmp_corr.pos[1]*param[1][1]+tmp_corr.pos[2]*param[2][1];
		m1.pos[2] = tmp_site.pos[2] + tmp_corr.pos[2]*param[2][2];
	}
	for(auto& m1 : B)
	{
		in>>tmp_site;
		getline(in,tmp);
		if(tmp_site.pos[0] - m1.pos[0] >= param[0][0]/2)
			tmp_corr.pos[0] = -1;
		else if(tmp_site.pos[0] - m1.pos[0] <= -param[0][0]/2)
			tmp_corr.pos[0] = 1;
		else
			tmp_corr.pos[0] = 0;
		if(tmp_site.pos[1] - m1.pos[1] >= param[1][1]/2)
			tmp_corr.pos[1] = -1;
		else if(tmp_site.pos[1] - m1.pos[1] <= -param[1][1]/2)
			tmp_corr.pos[1] = 1;
		else
			tmp_corr.pos[1] = 0;
		if(tmp_site.pos[2] - m1.pos[2] >= param[2][2]/2)
			tmp_corr.pos[2] = -1;
		else if(tmp_site.pos[2] - m1.pos[2] <= -param[2][2]/2)
			tmp_corr.pos[2] = 1;
		else
			tmp_corr.pos[2] = 0;
		m1.pos[0] = tmp_site.pos[0] + tmp_corr.pos[0]*param[0][0]+tmp_corr.pos[1]*param[1][0]+tmp_corr.pos[2]*param[2][0];
		m1.pos[1] = tmp_site.pos[1] + tmp_corr.pos[1]*param[1][1]+tmp_corr.pos[2]*param[2][1];
		m1.pos[2] = tmp_site.pos[2] + tmp_corr.pos[2]*param[2][2];
	}
	for(auto& m1 : C)
	{
		in>>tmp_site;
		getline(in,tmp);
		if(tmp_site.pos[0] - m1.pos[0] >= param[0][0]/2)
			tmp_corr.pos[0] = -1;
		else if(tmp_site.pos[0] - m1.pos[0] <= -param[0][0]/2)
			tmp_corr.pos[0] = 1;
		else
			tmp_corr.pos[0] = 0;
		if(tmp_site.pos[1] - m1.pos[1] >= param[1][1]/2)
			tmp_corr.pos[1] = -1;
		else if(tmp_site.pos[1] - m1.pos[1] <= -param[1][1]/2)
			tmp_corr.pos[1] = 1;
		else
			tmp_corr.pos[1] = 0;
		if(tmp_site.pos[2] - m1.pos[2] >= param[2][2]/2)
			tmp_corr.pos[2] = -1;
		else if(tmp_site.pos[2] - m1.pos[2] <= -param[2][2]/2)
			tmp_corr.pos[2] = 1;
		else
			tmp_corr.pos[2] = 0;
		m1.pos[0] = tmp_site.pos[0] + tmp_corr.pos[0]*param[0][0]+tmp_corr.pos[1]*param[1][0]+tmp_corr.pos[2]*param[2][0];
		m1.pos[1] = tmp_site.pos[1] + tmp_corr.pos[1]*param[1][1]+tmp_corr.pos[2]*param[2][1];
		m1.pos[2] = tmp_site.pos[2] + tmp_corr.pos[2]*param[2][2];
	}
}

void cell :: regist_oct()
{
	site tmp;
	double dis_aa = 4.7, dis_ab = 4.2, dis_bc = 2.8, dis;
	double xx2 = param[0][0]/2, yy2 = param[1][1]/2, zz2 = param[2][2]/2;
	// regist B site atom
	for(size_t t1=0; t1<num_b; t1++)
		oct[t1].B = B[t1];
	// regist A site atoms, use A-B < 4.2
	for(auto& m1 : oct)
	for(size_t t2=0; t2<num_a; t2++)
	{
		tmp.pos[0] = tmp.pos[1] = tmp.pos[2] = 0;
		if(m1.B.pos[0] - A[t2].pos[0] >= xx2)
			tmp.pos[0] = 1;
		else if(m1.B.pos[0] - A[t2].pos[0] <= -xx2)
			tmp.pos[0] = -1;
		if(m1.B.pos[1] - A[t2].pos[1] >= yy2)
			tmp.pos[1] = 1;
		else if(m1.B.pos[1] - A[t2].pos[1] <= -yy2)
			tmp.pos[1] = -1;
		if(m1.B.pos[2] - A[t2].pos[2] >= zz2)
			tmp.pos[2] = 1;
		else if(m1.B.pos[2] - A[t2].pos[2] <= -zz2)
			tmp.pos[2] = -1;
		dis = sqrt(pow((A[t2].pos[0]+tmp.pos[0]*param[0][0]+tmp.pos[1]*param[1][0]+tmp.pos[2]*param[2][0]-m1.B.pos[0]),2) + pow((A[t2].pos[1]+tmp.pos[1]*param[1][1]+tmp.pos[2]*param[2][1]-m1.B.pos[1]),2) + pow((A[t2].pos[2]+tmp.pos[2]*param[2][2]-m1.B.pos[2]),2));
		if(dis < dis_ab)
		{
			m1.ind_A.push_back(t2);
			m1.corr_A.push_back(tmp);
			m1.A.push_back(A[t2]);
			m1.A.back().pos[0] += tmp.pos[0]*param[0][0]+tmp.pos[1]*param[1][0]+tmp.pos[2]*param[2][0];
			m1.A.back().pos[1] += tmp.pos[1]*param[1][1]+tmp.pos[2]*param[2][1];
			m1.A.back().pos[2] += tmp.pos[2]*param[2][2];
		}
	}
	// regist C site atoms, use B-C < 2.8
	for(auto& m1 : oct)
	for(size_t t2=0; t2<num_c; t2++)
	{
		tmp.pos[0] = tmp.pos[1] = tmp.pos[2] = 0;
		if(m1.B.pos[0] - C[t2].pos[0] >= xx2)
			tmp.pos[0] = 1;
		else if(m1.B.pos[0] - C[t2].pos[0] <= -xx2)
			tmp.pos[0] = -1;
		if(m1.B.pos[1] - C[t2].pos[1] >= yy2)
			tmp.pos[1] = 1;
		else if(m1.B.pos[1] - C[t2].pos[1] <= -yy2)
			tmp.pos[1] = -1;
		if(m1.B.pos[2] - C[t2].pos[2] >= zz2)
			tmp.pos[2] = 1;
		else if(m1.B.pos[2] - C[t2].pos[2] <= -zz2)
			tmp.pos[2] = -1;
		dis = sqrt(pow((C[t2].pos[0]+tmp.pos[0]*param[0][0]+tmp.pos[1]*param[1][0]+tmp.pos[2]*param[2][0]-m1.B.pos[0]),2) + pow((C[t2].pos[1]+tmp.pos[1]*param[1][1]+tmp.pos[2]*param[2][1]-m1.B.pos[1]),2) + pow((C[t2].pos[2]+tmp.pos[2]*param[2][2]-m1.B.pos[2]),2));
		if(dis < dis_bc)
		{
			m1.ind_C.push_back(t2);
			m1.corr_C.push_back(tmp);
			m1.C.push_back(C[t2]);
			m1.C.back().pos[0] += tmp.pos[0]*param[0][0]+tmp.pos[1]*param[1][0]+tmp.pos[2]*param[2][0];
			m1.C.back().pos[1] += tmp.pos[1]*param[1][1]+tmp.pos[2]*param[2][1];
			m1.C.back().pos[2] += tmp.pos[2]*param[2][2];
		}
	}
//	cout<<m1.C.size()<<endl;
/*
	for(auto& m1 : oct)
	{
		cout<<"Fe "<<m1.B<<endl;
		for(auto& m2 : m1.C)
			cout<<"O "<<m2<<endl;
	}
*/
}

void cell :: rebuild_oct()
{
	for(size_t t1=0; t1<num_b; t1++)
	{
		oct[t1].B = B[t1];
		for(size_t t2=0; t2<oct[t1].A.size(); t2++)
		{
			oct[t1].A[t2].pos[0] = A[oct[t1].ind_A[t2]].pos[0] + oct[t1].corr_A[t2].pos[0]*param[0][0]+oct[t1].corr_A[t2].pos[1]*param[1][0]+oct[t1].corr_A[t2].pos[2]*param[2][0];
			oct[t1].A[t2].pos[1] = A[oct[t1].ind_A[t2]].pos[1] + oct[t1].corr_A[t2].pos[1]*param[1][1]+oct[t1].corr_A[t2].pos[2]*param[2][1];
			oct[t1].A[t2].pos[2] = A[oct[t1].ind_A[t2]].pos[2] + oct[t1].corr_A[t2].pos[2]*param[2][2];
		}
		for(size_t t2=0; t2<oct[t1].C.size(); t2++)
		{
			oct[t1].C[t2].pos[0] = C[oct[t1].ind_C[t2]].pos[0] + oct[t1].corr_C[t2].pos[0]*param[0][0]+oct[t1].corr_C[t2].pos[1]*param[1][0]+oct[t1].corr_C[t2].pos[2]*param[2][0];
			oct[t1].C[t2].pos[1] = C[oct[t1].ind_C[t2]].pos[1] + oct[t1].corr_C[t2].pos[1]*param[1][1]+oct[t1].corr_C[t2].pos[2]*param[2][1];
			oct[t1].C[t2].pos[2] = C[oct[t1].ind_C[t2]].pos[2] + oct[t1].corr_C[t2].pos[2]*param[2][2];
		}
	}

	for(auto& m1 : oct)
	{
		cout<<"Fe "<<m1.B<<endl;
		for(auto& m2 : m1.A)
			cout<<"Bi "<<m2<<endl;
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
