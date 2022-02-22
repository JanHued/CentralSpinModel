/*
 * spindyn_sec.cpp

 */
#include <iostream>
#include <fstream>
#define _USE_MATH_DEFINES
#include <math.h>
#include <vector>
#include <chrono>
#include <random>
#include <iomanip>

using namespace std;
//#include <iomanip>
#include <string>

unsigned seed1 = chrono::system_clock::now().time_since_epoch().count();
	//default_random_engine generator;
	//random_device rd; //Seed fuer Mersenne-Twister
	mt19937 gen(seed1); //Mersenne-Twister
	//uniform_real_distribution<double> gleich(0, 1); //Gleichverteilung
	normal_distribution<double> normal(0.0, 1./2.); //Gaussverteilung

vector<double> product(vector<double> a, vector<double> b){
	vector<double> result={0,0,0};
	result[0]=a[1]*b[2]-a[2]*b[1];
	result[1]=a[2]*b[0]-a[0]*b[2];
	result[2]=a[0]*b[1]-a[1]*b[0];
	return result;
}

vector<double> operator * (const double& faktor, const vector<double>& v1)
{
	vector<double> result=v1;
	for (unsigned int i = 0; i < result.size(); ++i)
	{
		result[i] *= faktor;
	}
	return result;
}

void operator *= (const double& faktor, vector<double>& v1)
{
	for (unsigned int i = 0; i < v1.size(); ++i)
	{
		v1[i] *= faktor;
	}
}

vector<double> operator + (const vector<double>& v1, const vector<double>& v2)
{
	vector<double> result=v1;
	for (unsigned int i = 0; i < result.size(); ++i)
	{
		result[i] += v2[i];
	}
	return result;
}

vector<double> operator - (const vector<double>& v1, const vector<double>& v2)
{
	vector<double> result=v1;
	for (unsigned int i = 0; i < result.size(); ++i)
	{
		result[i] -= v2[i];
	}
	return result;
}

void operator += (vector<double>& v1, const vector<double>& v2)
{
	for (unsigned int i = 0; i < v1.size(); ++i)
	{
		v1[i] += v2[i];
	}
}

double operator * (const vector<double>& v1, const vector<double>& v2)
{
	double result=0;
	for(unsigned int i=0;i<3;i++){
		result+=v1[i]*v2[i];
	}
	return result;
}

void out(vector<double> v){
	int size=v.size();
	for(int i=0;i<size;i++){
	cout << v[i] << endl;
	}
	cout << endl;
}

void precession(vector<double> &res, const vector<double>& S0,const vector<double>& n, double omega, double t){
	vector<double> tmp=(S0*n)*n;
	vector<double> tmp2=S0-tmp;
	res=tmp+cos(omega*t)*(tmp2)-sin(omega*t)*product(tmp2,n);
}

vector<double> precession(const vector<double>& S0,const vector<double>& n, double omega, double t){
	vector<double>res=(S0*n)*n;
	double arg=omega*t;
	vector<double>res1=S0+(-1)*res;
	res+=cos(arg)*res1-sin(arg)*product(res1,n);
	return res;
}

void couplings_equidist(vector<double>& J,int N, double Jq){
	for(int i=0;i<N;i++){
			J[i]=sqrt(6.*N/(2.*N*N+3.*N+1.))*(N-i)/N*Jq;
		}
}

void couplings_exp(vector<double>& J, int N, double x){
	double norm;
	norm=(exp(-2*x*(1+1./N))-1)/(exp(-2.*x/N)-1)-1.0;
	norm=sqrt(norm);
	for(int i=1;i<=N;i++){
		J[i]=1./norm*exp(-double(i)*x/double(N));
	}
}

void couplings_exp_gamma(vector<double>& J, int N, double gamma){
	double A0=sqrt((1-exp(-2*gamma))/(1-exp(-2*gamma*N)));
	for(int i=0;i<N;i++){
		J[i]=A0*exp(-double(i)*gamma);
	}
}

void dynamic_fkt(int N,int max_configs, double h, int max_step, string name, double x, double gamma, int coupling){
	//Central-Spin
	vector<double> S(3,0.);
	vector<double> S0(3,0.);
	//Mean
	vector<double> mean(max_step,0.);
	//Bath-Spins
	vector<vector<double> > B(N,vector<double>(3,0.));
	vector<vector<double> > Bnext(N,vector<double>(3,0.));
	vector<vector<double> > B0(N,vector<double>(3,0.));
	//couplings
	vector<double> J(N,0.);
	//overhauser
	vector<double> A(3,0.);
	double omega_CS=0;
	vector<double> omega_bath(N,0.);
	vector<double> n_A(3,0.);
	vector<double> n_CS(3,0.);
	double A_abs;
	double S_abs;
	vector<double> start;
	double datum=0;
	int ctr=0;

	fstream startwerte("badspins_start.txt",ios::in);
	while(startwerte>>datum){
		start.push_back(datum);
	}
	startwerte.close();

	//Test
	double energy1=0;
	double energy2=0;

	//Initialize couplings
	switch(coupling){
	case 1: couplings_equidist(J,N,1.);break;
	case 2: couplings_exp(J,N,x);break;
	case 3: couplings_exp_gamma(J,N,gamma);break;
	}

	//Calculate Dynamic
	for(int config=0;config<max_configs;config++){
//		if(config%(max_configs/10)==0){
//			cout << "Fortschritt Funktion direkt: " << int(config*10/max_configs) << "/10" << endl;
//		}
		//Initialize random Config
//		for(int n=0;n<N;n++){
//			for(int comp=0;comp<3;comp++){
//				B[n][comp]=normal(gen);
//			}
//		}
		ctr=0;
		//Init with known values
		for(int n=0;n<N;n++){
			for(int comp=0;comp<3;comp++){
				B[n][comp]=start[ctr];
				ctr++;
			}
		}

		//rnd CS
//		S[0]=normal(gen);
//		S[1]=normal(gen);
//		S[2]=normal(gen);
		S[0]=0.;
		S[1]=0.;
		S[2]=1.;

		//values at t=0
		S0=S;
		B0=B;

		//half step of overhauser
//		S_abs=sqrt(S*S);
//		n_CS=1./S_abs*S;
//		omega_bath=S_abs*J;
//		for(int i=0;i<N;i++){
//			precession(B[i],B[i],n_CS,omega_bath[i],h/2.);
//		}

		for(int step=0;step<max_step;step++){
			//Test
//			if(step==0){
//				for(int i=0;i<N;i++){
//					energy1+=J[i]*S*B[i];
//				}
//			}
//			if(step==max_step-1){
//				for(int i=0;i<N;i++){
//					energy2+=J[i]*S*B[i];
//				}
//				cout << "S(t): E(0)=" << energy1 << endl;
//				cout << "S(t): dE mit h=" << h << ": " << abs((energy2-energy1)/energy1) << endl;
//			}

			//ensemble mean
//			mean[step]+=S[2]*S0[2]/max_configs;
			A.assign(3,0.);
			//Calculate Overhauser-Field
			for(int i=0;i<N;i++){
				A=A+J[i]*B[i];
			}
			A_abs=sqrt(A*A);
			S_abs=sqrt(S*S);
			n_A=1./A_abs*A;
			n_CS=1./S_abs*S;
			omega_CS=A_abs;
			omega_bath=S_abs*J;
			//calc precession of CS
			precession(S,S,n_A,omega_CS,h);
			//calc precession of bathspins
			for(int i=0;i<N;i++){
				precession(B[i],B[i],n_CS,omega_bath[i],h);
			}
			if(step==max_step-1){
				cout << "Merk. mit h=" << h << ":   Sz(100) = " << setprecision(10) << S[2] << endl;
			}
		}
	}
//	fstream daten(name,ios::out);
//	for(int step=0;step<max_step;step++){
//		daten << h*step << "\t" << mean[step] << endl;
//	}
}

int main(int argN, char** args){
	double tmax=100;
	int N=100;
	int ensembles=10e5;
	double h=0.01;
	double x=1.;
	double gamma=0.01;
	string zusatz = "test";
	int coupling = 0;			// 1: equi		2: exp		3: gamma

	if(argN == 8)
	{
		ensembles=atoi(args[1]);
		x=atoi(args[2]);
		h=atof(args[3]);
		N=atoi(args[4]);
		zusatz=args[5];
		coupling=atoi(args[6]);
		gamma=atof(args[7]);
	}
	else{
		cout << "Falsche Anzahl Parameter. Verlangt: 7" << endl;
		return -1;
	}
	int steps=tmax/h;

	string name="conserv_test";

	//string name="spindyn_S_";
//	string name="testSequi_";
//	name += zusatz;
	name += ".txt";

	dynamic_fkt(N,ensembles,h,steps,name, x, gamma, coupling);

	return 0;
}
