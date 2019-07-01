#include <cstring>
#include <sstream>
#include <iomanip>
#include <cmath>
#include <fstream>
#include "cell.h"

#ifdef __SPECTRA__
#include <fftw3.h>
#endif

using namespace std;

void cell :: init(ifstream& in)
{
	string label_iter = "ITEM: TIMESTEP";
	string label_num = "ITEM: NUMBER OF ATOMS";
	string tmp;
	// get number of iteration
	tot_step = 0;
	while(!in.eof())
	{
		getline(in,tmp);
		if(tmp.find(label_iter) != string::npos)
			tot_step++;
	}
	in.clear(); in.seekg(ios::beg);
	polarization.resize(tot_step);
	polarization_stderr.resize(tot_step);
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
	string label_cell = "xlo xhi";
	string label_tri = "xy xz yz";
	string label_atom = "Atoms";
	stringstream ss;
	string tmp;
	// find cell parameter
	getline(in,tmp);
	while(tmp.find(label_cell) == string::npos)
		getline(in,tmp);
	ss<<(tmp);
	ss>>shift[0]>>param[0][0];
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
	// check if triclinic or not
	getline(in,tmp);
	if(tmp.find(label_tri) != string::npos)
	{
		ss.str(""); ss.clear();
		ss<<(tmp);
		ss>>param[1][0]>>param[2][0]>>param[2][1];
	}
	// get position of atoms
	getline(in,tmp);
	while(tmp.find(label_atom) == string::npos)
		getline(in,tmp);
	getline(in,tmp);
	for(auto& m1 : A)
	{
		in>>tmp>>tmp>>tmp>>tmp>>m1;
		getline(in,tmp);
	}
	for(auto& m1 : B)
	{
		in>>tmp>>tmp>>tmp>>tmp>>m1;
		getline(in,tmp);
	}
	for(auto& m1 : C)
	{
		in>>tmp>>tmp>>tmp>>tmp>>m1;
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
		// adjust based on https://lammps.sandia.gov/doc/Howto_triclinic.html
		if (param[2][1]>0) param[1][1] -= param[2][1];
		if (param[2][1]<0) shift[1]    -= param[2][1];
		if (max(max(param[1][0],param[2][0]),param[1][0]+param[2][0])>0) param[0][0] -= max(max(param[1][0],param[2][0]),param[1][0]+param[2][0]);
		if (min(min(param[1][0],param[2][0]),param[1][0]+param[2][0])<0) shift[0]    -= min(min(param[1][0],param[2][0]),param[1][0]+param[2][0]);
		// end adjust
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
	ofstream out;
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
	// save neighbor to file
	out.open("neighbor_list.dat");
	for(size_t t1=0; t1<num_b; t1++)
	{
		out<<t1;
		out<<' '<<oct[t1].ind_A.size()<<' '<<oct[t1].ind_C.size();
		for(auto& m2 : oct[t1].ind_A)
			out<<' '<<m2;
		for(auto& m2 : oct[t1].ind_C)
			out<<' '<<m2;
		out<<endl;
	}
	out.close();
/*
	for(auto& m1 : oct)
	{
		cout<<"Fe "<<m1.B<<endl;
		for(auto& m2 : m1.C)
			cout<<"O "<<m2<<endl;
	}
*/
}

int cell :: recover_oct()
{
	ifstream in;
	in.open("neighbor_list.dat");
	double xx2 = param[0][0]/2, yy2 = param[1][1]/2, zz2 = param[2][2]/2;
	int nn_a, nn_c;
	string tmp;
	if(in.good())
	{
		// get index
		for(size_t t1=0; t1<num_b; t1++)
		{
			// site B
			oct[t1].B = B[t1];
			in>>tmp>>nn_a>>nn_c;
			oct[t1].A.resize(nn_a);
			oct[t1].corr_A.resize(nn_a);
			oct[t1].ind_A.resize(nn_a);
			oct[t1].C.resize(nn_c);
			oct[t1].corr_C.resize(nn_c);
			oct[t1].ind_C.resize(nn_c);
			// site A
			for(size_t t2=0; t2<nn_a; t2++)
				in>>oct[t1].ind_A[t2];
			// site C
			for(size_t t2=0; t2<nn_c; t2++)
				in>>oct[t1].ind_C[t2];
			getline(in,tmp);
		}
		// get correction
		for(auto& m1 : oct)
		{
			// site A
			for(size_t t2=0; t2<m1.ind_A.size(); t2++)
			{
				m1.corr_A[t2].clear();
				if(m1.B.pos[0] - A[m1.ind_A[t2]].pos[0] >=  xx2) m1.corr_A[t2].pos[0] =  1;
				else if(m1.B.pos[0] - A[m1.ind_A[t2]].pos[0] <= -xx2) m1.corr_A[t2].pos[0] = -1;
				if(m1.B.pos[1] - A[m1.ind_A[t2]].pos[1] >=  yy2) m1.corr_A[t2].pos[1] =  1;
				else if(m1.B.pos[1] - A[m1.ind_A[t2]].pos[1] <= -yy2) m1.corr_A[t2].pos[1] = -1;
				if(m1.B.pos[2] - A[m1.ind_A[t2]].pos[2] >=  zz2) m1.corr_A[t2].pos[2] =  1;
				else if(m1.B.pos[2] - A[m1.ind_A[t2]].pos[2] <= -zz2) m1.corr_A[t2].pos[2] = -1;
			}
			// site C
			for(size_t t2=0; t2<m1.ind_C.size(); t2++)
			{
				m1.corr_C[t2].clear();
				if(m1.B.pos[0] - C[m1.ind_C[t2]].pos[0] >=  xx2) m1.corr_C[t2].pos[0] =  1;
				else if(m1.B.pos[0] - C[m1.ind_C[t2]].pos[0] <= -xx2) m1.corr_C[t2].pos[0] = -1;
				if(m1.B.pos[1] - C[m1.ind_C[t2]].pos[1] >=  yy2) m1.corr_C[t2].pos[1] =  1;
				else if(m1.B.pos[1] - C[m1.ind_C[t2]].pos[1] <= -yy2) m1.corr_C[t2].pos[1] = -1;
				if(m1.B.pos[2] - C[m1.ind_C[t2]].pos[2] >=  zz2) m1.corr_C[t2].pos[2] =  1;
				else if(m1.B.pos[2] - C[m1.ind_C[t2]].pos[2] <= -zz2) m1.corr_C[t2].pos[2] = -1;
			}
		}
		return 1;
	}
	else
		return 0;
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
/*
	for(auto& m1 : oct)
	{
		cout<<"Fe "<<m1.B<<endl;
		for(auto& m2 : m1.A)
			cout<<"Bi "<<m2<<endl;
	}
*/
}

void cell :: ave_p(int iter)
{
	site res, tmp;
	res.clear();
	polarization[iter].clear();
	polarization_stderr[iter].clear();
	for(auto& m1 : oct)
	{
		res.clear();
		// add B
		res = m1.B*q_b;
		//cout<<"Fe "<<m1.B<<endl;
		// add A
		tmp.clear();
		for(auto& m2 : m1.A)
		{
			tmp = tmp + m2;
			//cout<<"Bi "<<m2<<endl;
		}
		tmp = tmp*q_a/m1.A.size();
		res = res + tmp;
		// add C
		tmp.clear();
		for(auto& m2 : m1.C)
		{
			tmp = tmp + m2;
			//cout<<"O  "<<m2<<endl;
		}
		tmp = tmp*q_c*3/m1.C.size();
		res = res + tmp;
		// add to output
		polarization[iter] = polarization[iter] + res;
		polarization_stderr[iter].pos[0] += res.pos[0]*res.pos[0];
		polarization_stderr[iter].pos[1] += res.pos[1]*res.pos[1];
		polarization_stderr[iter].pos[2] += res.pos[2]*res.pos[2];
		//cout<<endl;
	}
	polarization[iter] = polarization[iter] / num_b;
	polarization_stderr[iter].pos[0] = sqrt(polarization_stderr[iter].pos[0]/num_b - polarization[iter].pos[0]*polarization[iter].pos[0]);
	polarization_stderr[iter].pos[1] = sqrt(polarization_stderr[iter].pos[1]/num_b - polarization[iter].pos[1]*polarization[iter].pos[1]);
	polarization_stderr[iter].pos[2] = sqrt(polarization_stderr[iter].pos[2]/num_b - polarization[iter].pos[2]*polarization[iter].pos[2]);
}

#ifdef __SPECTRA__
void cell :: init_spectra()
{
	dw = 1.0/(dt*(tot_step-1));
	freq.resize(tot_step);
	p_w.resize(tot_step);
	spectra.resize(tot_step);
	traj = new double*[tot_step];
	for(size_t t1=0; t1<tot_step; t1++)
		traj[t1] = new double[(num_a+num_b+num_c)*3];
}

void cell:: save_traj(int iter)
{
	freq[iter] = iter*dw;
	for(size_t t1=0; t1<num_a; t1++)
	{
		traj[iter][t1*3+0] = A[t1].pos[0];
		traj[iter][t1*3+1] = A[t1].pos[1];
		traj[iter][t1*3+2] = A[t1].pos[2];
	}
	for(size_t t1=0; t1<num_b; t1++)
	{
		traj[iter][(num_a+t1)*3+0] = B[t1].pos[0];
		traj[iter][(num_a+t1)*3+1] = B[t1].pos[1];
		traj[iter][(num_a+t1)*3+2] = B[t1].pos[2];
	}
	for(size_t t1=0; t1<num_c; t1++)
	{
		traj[iter][(num_a+num_b+t1)*3+0] = C[t1].pos[0];
		traj[iter][(num_a+num_b+t1)*3+1] = C[t1].pos[1];
		traj[iter][(num_a+num_b+t1)*3+2] = C[t1].pos[2];
	}
}

void cell :: get_p_w()
{
	double mean;
	// initialize fft
	fftw_complex *input, *output;
	fftw_plan p0;

	input  = new fftw_complex[tot_step];
	output = new fftw_complex[tot_step];
	p0 = fftw_plan_dft_1d(tot_step, input, output, -1, FFTW_ESTIMATE);

	// calculate p_w
	for(size_t t1=0; t1<3; t1++)
	{
		mean = 0;
		for(size_t t2=0; t2<tot_step; t2++)
			mean += polarization[t2].pos[t1];
		mean /= tot_step;
		for(size_t t2=0; t2<tot_step; t2++)
		{
			input[t2][0] = polarization[t2].pos[t1] - mean;
			input[t2][1] = 0;
		}
		fftw_execute(p0);
		for(size_t t2=0; t2<tot_step; t2++)
			p_w[t2].pos[t1] = sqrt(output[t2][0]*output[t2][0]+output[t2][1]*output[t2][1]);
	}
}

void cell :: get_spectra()
{
	double mean;
	// initialize fft
	fftw_complex *input, *output;
	fftw_plan p0;

	input  = new fftw_complex[tot_step];
	output = new fftw_complex[tot_step];
	p0 = fftw_plan_dft_1d(tot_step, input, output, -1, FFTW_ESTIMATE);

	for(size_t t1=0; t1<tot_step; t1++)
		spectra[t1] = 0;
	// loop for each degree of freedom(DOF)
	for(size_t t1=0; t1<(num_a+num_b+num_c)*3; t1++)
	{
		mean = 0;
		for(size_t t2=0; t2<tot_step; t2++)
			mean += traj[t2][t1];
		mean /= tot_step;
		for(size_t t2=0; t2<tot_step; t2++)
		{
			input[t2][0] = traj[t2][t1] - mean;
			input[t2][1] = 0;
		}
		fftw_execute(p0);
		for(size_t t2=0; t2<tot_step; t2++)
			spectra[t2] += (output[t2][0]*output[t2][0]+output[t2][1]*output[t2][1]);
	}
	for(size_t t1=0; t1<tot_step; t1++)
		spectra[t1] *= freq[t1]*freq[t1]/3/(num_a+num_b+num_c);
}
#endif



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
