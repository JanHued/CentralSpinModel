/*
 * continous_fraction.cpp
 *
 *  Created on: 15.02.2016
 *      Author: huedepohl
 */

#define DEBUG 0
#include <iostream>
#include <fstream>
#define _USE_MATH_DEFINES
//#include <math.h>
#include <vector>
#include <iomanip>
#include <string>
//#include <eigen3/Eigen/Dense>
//#include <functional>
#include <cassert>
#include <sstream>
#include <complex>
#include <cmath>
using namespace std;
//using namespace Eigen;
//using Eigen::MatrixXd;

constexpr complex<double> I{0,1};

void out(const vector<double>& v){
	int size=v.size();
	for(int i=0;i<size;i++){
		cout << v[i] << endl;
	}
}

void out(const vector<vector<double>>& v){
	int a=v.size();
	int b=v[0].size();
	for(int i=0;i<b;i++){
		for(int j=0;j<a;j++){
			cout << v[j][i] << "\t";
		}
		cout << endl;
	}
}

void couplings_equidist(vector<double>& J){
	int N=J.size();
	for(int i=0;i<N;i++){
			J[i]=sqrt(6.*N/(2.*N*N+3.*N+1.))*(N-i)/N;
		}
}

void couplings_gamma(vector<double>& J, const double gamma){
	assert(gamma<1);
	const int N=J.size();
	assert(N>1);
	const double A0=sqrt((1.-exp(-2.*gamma))/(1.-exp(-2.*gamma*N)));
	for(int i=0;i<N;i++){
		J[i]=A0*exp(-double(i)*gamma);
	}
}

void read(const string alpha,const string beta,vector<vector<double>>&coeff,const int depth){
	long double x=0.;
	fstream coeffa(alpha,ios::in);
	fstream coeffb(beta,ios::in);
	while(coeff[0].size()<depth){
		coeffa>>x;
		coeff[0].push_back(x);
		x=0.;
		coeffb>>x;
		coeff[1].push_back(x);
	}
}

void continued_fraction(const vector<vector<double>>&coeff,const complex<double>w,complex<double>& res){
	int size=coeff[0].size();
	res=w-coeff[0][size-1]-coeff[1][size-1];
	for(unsigned i=1;i<size;++i){
		res=w-coeff[0][size-1-i]-coeff[1][size-1-i]/res;
	}
	res=1./res;
}

void resolvente(const vector<vector<double>>&coeff,vector<vector<double>>&res,const double w_min,const double w_max,const int steps, const int depth,double delta){
	complex<double>R=0.;
	double invpi=-1./M_PI;
	double dw=(w_max-w_min)/steps;
	complex<double>w;
	w=w_min+delta*I;
	for(unsigned i=0;i<steps;++i){
		w+=dw;
		continued_fraction(coeff,w,R);
		res[0].push_back(w.real());
		res[1].push_back(invpi*R.imag());
	}
}


int main(int argN,char**args) {
	assert(argN==7);

	int N=atoi(args[1]);
	char* spins=args[1];
	int couptype=atoi(args[2]);	// 1: equi	2: gamma
	char* g=args[3];
	double gamma=atof(args[3]);
	int depth=atoi(args[4]);
	char* depth_=args[4];
	int steps=atoi(args[5]);
	double delta=atof(args[6]);

	vector<vector<double>> coeff(2,vector<double>());
	vector<vector<double>> res(2,vector<double>());

	string filealpha;
	string filebeta;
	string datenfile;
	stringstream sdoub;
	stringstream sdelta;
	sdoub << setprecision(2) << g;
	sdelta << delta;
	
	vector<double> couplings(N,0.);
	switch(couptype){
		case 1: couplings_equidist(couplings);
				filealpha="equi/alpha_N"+string(spins)+".txt";
				filebeta="equi/beta_N"+string(spins)+".txt";
				datenfile="density_equi_N"+string(spins)+"_depth_"+string(depth_)+"delta_"+sdelta.str()+".txt";
				break;
		case 2: if(N==1){
					filealpha="gamma/alpha_gamma"+sdoub.str()+"_Ninfty.txt";
					filebeta="gamma/beta_gamma"+sdoub.str()+"_Ninfty.txt";
					datenfile="density_gamma"+sdoub.str()+"_Ninfty_depth_"+string(depth_)+"delta_"+sdelta.str()+".txt";
				}
				else{
					couplings_gamma(couplings,gamma);
					filealpha="gamma/alpha_gamma"+sdoub.str()+"_N"+string(spins)+".txt";
					filebeta="gamma/beta_gamma"+sdoub.str()+"_N"+string(spins)+".txt";
					datenfile="density_gamma"+sdoub.str()+"_N"+string(spins)+"_depth_"+string(depth_)+"delta_"+sdelta.str()+".txt";
				}
					break;
	}
	

	//read coefficients alpha,beta
	read(filealpha,filebeta,coeff,depth);
	assert(coeff[0].size()==coeff[1].size());
	assert(depth<=coeff[0].size());

	double w_min=couplings[N-1];
	double w_max=couplings[0];
	if(N==1){
			w_min=0.;
			w_max=sqrt(1-exp(-2*gamma));
	}
	double w_tol=(w_max-w_min)/4;

	w_min-=w_tol;
	w_max+=w_tol;

	//only need beta^2
	for(unsigned i=0;i<coeff[0].size();++i){
		coeff[1][i]=coeff[1][i]*coeff[1][i];
	}
	resolvente(coeff,res,w_min,w_max,steps,depth,delta);
	fstream daten(datenfile,ios::out);
	for(unsigned i=0;i<res[0].size();++i){
		daten << res[0][i] << "\t" << res[1][i] << endl;
	}

	return 0;
}

