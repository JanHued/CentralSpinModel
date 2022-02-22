/*
 * rk_new.cpp
 *
 *  Created on: 02.02.2016
 *      Author: huedepohl
 */

#include <iostream>
#include <fstream>
#define _USE_MATH_DEFINES
#include <math.h>
#include <chrono>
#include <random>
#include <iomanip>
#include <functional>
using namespace std;
//#include <iomanip>
#include <string>
#include <cassert>
#include <gsl/gsl_odeiv2.h>
#include <gsl/gsl_errno.h>


unsigned seed1 = chrono::system_clock::now().time_since_epoch().count();
	//default_random_engine generator;
	//random_device rd; //Seed fuer Mersenne-Twister
	mt19937 gen(seed1); //Mersenne-Twister
	//uniform_real_distribution<double> gleich(0, 1); //Gleichverteilung
	normal_distribution<double> normal(0.0, 1./2.); //Gaussverteilung

struct parameter{
		int size; //Dim of ODEs
		double* couplings; //Couplings
		double B0; //B-Field value
		int B_function; //B-field function
	};

double Bfield(const double t,const double B0, const int function){
	assert(function>0);
	switch(function){
		case 1: return  B0; //B=const
	}
}

int rightside(double t, const double* y0, double* res, void* params){
	//y0[0],[1],[2] is central spin, rest bath spins

	parameter* p=(parameter *)params;
	double* J=p->couplings;
	int y0_size=p->size;
	double H=p->B0;
	int Bfunction=p->B_function;
	H=Bfield(t,H,Bfunction);

	int N=y0_size/3.-1;

	double A[3]={0.};
	//calculate overhauser field
	for(unsigned i=1;i<=N;++i){
		A[0]+=J[i-1]*y0[3*i];
		A[1]+=J[i-1]*y0[3*i+1];
		A[2]+=J[i-1]*y0[3*i+2];
	}
	//right side of CS with B-Feld H in z-dir.
	res[0]=A[1]*y0[2]-A[2]*y0[1]+H*y0[1];
	res[1]=A[2]*y0[0]-A[0]*y0[2]-H*y0[0];
	res[2]=A[0]*y0[1]-A[1]*y0[0];

	//right side of BS
	for(unsigned i=1;i<=N;++i){
		res[3*i]=J[i-1]*(y0[1]*y0[3*i+2]-y0[2]*y0[3*i+1]);
		res[3*i+1]=J[i-1]*(y0[2]*y0[3*i]-y0[0]*y0[3*i+2]);
		res[3*i+2]=J[i-1]*(y0[0]*y0[3*i+1]-y0[1]*y0[3*i]);
	}
	return GSL_SUCCESS;
}

void out(const double* v, const int v_size){
	for(int i=0;i<v_size;i++){
	cout << v[i] << endl;
	}
	cout << endl;
}


void couplings_equidist(double* J, const int N){
	double Jq=1.;
	for(int i=0;i<N;i++){
			J[i]=sqrt(6.*N/(2.*N*N+3.*N+1.))*(N-i)/N*Jq;
		}
}

void couplings_exp_gamma(double* J, const int N, const double gamma){
	double A0=sqrt((1-exp(-2*gamma))/(1-exp(-2*gamma*N)));
	for(int i=0;i<N;i++){
		J[i]=A0*exp(-double(i)*gamma);
	}
}

void couplings_lin(double* J, const int N){
	double J0=sqrt(1./(1./6.*N*(N+1)*(2*N+1)));
	for(int i=0;i<N;i++){
		J[i]=J0*(i+1);
	}
}

void couplings_posneg(double* J, const int J_size, const double gamma){
	assert(gamma<1);
	const int N=J_size;
	assert(N%2==0);
	assert(N>1);
	const double A0=sqrt((1.-exp(-2.*gamma))/(1.-exp(-gamma*N)));
	for(int i=0;i<N/2;i++){
		J[i]=A0/sqrt(2)*exp(-double(i)*gamma);
	}
	for(int i=0;i<N/2;i++){
		J[i+N/2]=-A0/sqrt(2)*exp(-double(i)*gamma);
	}
}

void flip_spin(double* y0, const int dir){
	// flip spin in direction dir  --> pulse
	double spinval=sqrt(y0[0]*y0[0]+y0[1]*y0[1]+y0[2]*y0[2]);
	y0[0]=0.;
	y0[1]=0.;
	y0[2]=0.;
	y0[dir]=spinval;
}

void pi_pulse(double* y0){
	y0[0]=-y0[0];
	y0[2]-y0[2];
}

void abs_pulse(double* y0, const int dir){
	y0[dir]=abs(y0[dir]);
}


void dynamic_RK(const int N, const int max_configs, const int tmax, const string name, const double gamma, const int coupling, const double H, const int datapoints, const int B_function, const double t_pulse){

	const double delta = double(tmax)/datapoints;

	//Mean
	double* sz_mean=new double[datapoints];
	double* sx_mean=new double[datapoints];
	double* sy_mean=new double[datapoints];

	double* szsz_mean=new double[datapoints];
	double* sxsx_mean=new double[datapoints];
	double* sysy_mean=new double[datapoints];

	//steps befor pulse
	int pulsestep=int(t_pulse/delta);


	//data
	double* szsz_err=new double[datapoints];
	double szsz=0.;
	double* sxsx_err=new double[datapoints];
	double sxsx=0.;
	double* sysy_err=new double[datapoints];
	double sysy=0.;
	double sxsx_oldmean,sysy_oldmean,szsz_oldmean;

	//for envelope
	double* sysx0=new double[datapoints];

	for(unsigned i=0;i<datapoints;++i){
		sx_mean[i]=0.;
		sxsx_mean[i]=0.;
		sxsx_err[i]=0.;
		sy_mean[i]=0.;
		sysy_mean[i]=0.;
		sysy_err[i]=0.;
		sz_mean[i]=0.;
		szsz_mean[i]=0.;
		szsz_err[i]=0.;
		sysx0[i]=0.;
	}

	//initial condition of CS
	double S0[3];

	//couplings
	double* J=new double[N];

	//dimension of system (3 comp. for CS and 3 comp. for every BS)
	unsigned dim=3*N+3;

	double* y0=new double[dim];

	//Initialize couplings
	switch(coupling){
	case 1: couplings_equidist(J,N);break;
	case 2: couplings_exp_gamma(J,N,gamma);break;
	case 3: couplings_posneg(J,N,gamma);break;
	}
	fstream daten(name,ios::out);

	parameter params;
	params.size=dim;
	params.couplings=J;
	params.B0=H;
	params.B_function=B_function;
	parameter* p=&params;


	//initialize odes
	gsl_odeiv2_system ode = {rightside, NULL, dim, p};

	//rk4 = classic 4th order RK, rkf45 = explicit embedded Runge-Kutta-Fehlberg (4, 5) method.
	gsl_odeiv2_driver* d = gsl_odeiv2_driver_alloc_y_new(&ode,gsl_odeiv2_step_rkf45,1e-6,1e-6,0);

//	gsl_odeiv2_driver_set_hmin(d,1e-6);
//	gsl_odeiv2_driver_set_hmax(d,delta);
//	gsl_odeiv2_driver* d = gsl_odeiv2_driver_alloc_y_new(&ode,gsl_odeiv2_step_rk4,0.002,1e-6,0.)

	//Calculate Dynamic
	for(int config=1;config<=max_configs;config++){
		double t = 0;
		double ti=0;

		if(config%(max_configs/10)==0){
			cout << "Fortschritt alle DGLs: " << int(config*10/max_configs) << "/10" << endl;
		}

//		Initialize random Config
		for(unsigned i=0;i<dim;i++){
			y0[i]=normal(gen);
		}

//		flip_spin(y0,0);

		S0[0]=y0[0];
		S0[1]=y0[1];
		S0[2]=y0[2];

		for(int step=0;step<datapoints;step++){
			ti=(step+1)*delta;

//			if(pulsestep){
//				if(step%pulsestep==0&&step!=0){
//					flip_spin(y0,0);
//				}
//			}

			sxsx=y0[0]*S0[0];
			sysy=y0[1]*S0[1];
			szsz=y0[2]*S0[2];

			sx_mean[step]+=y0[0]/max_configs;
			sy_mean[step]+=y0[1]/max_configs;
			sz_mean[step]+=y0[2]/max_configs;

			sxsx_mean[step]+=sxsx/max_configs;
			sysy_mean[step]+=sysy/max_configs;
			szsz_mean[step]+=szsz/max_configs;

			sysx0[step]+=(y0[1]*S0[0])/max_configs;


//			//averaging S_0 with variance (Welford-Alg)
//			sxsx_oldmean=sxsx_mean[step];
//			sxsx_mean[step]+=(sxsx-sxsx_mean[step])/config;
//			sxsx_err[step]+=(sxsx-sxsx_mean[step])*(sxsx-sxsx_oldmean);
//
//			sysy_oldmean=sysy_mean[step];
//			sysy_mean[step]+=(sysy-sysy_mean[step])/config;
//			sysy_err[step]+=(sysy-sysy_mean[step])*(sysy-sysy_oldmean);
//
//			szsz_oldmean=szsz_mean[step];
//			szsz_mean[step]+=(szsz-szsz_mean[step])/config;
//			szsz_err[step]+=(szsz-szsz_mean[step])*(szsz-szsz_oldmean);

            //time evolution
            int status = gsl_odeiv2_driver_apply(d,&t,ti,y0);

            if(status != GSL_SUCCESS) {
                cerr << "GSL Error: return value = " << status << endl;
                exit(1);
            }
		}
	}
    gsl_odeiv2_driver_free(d);

// 	save data of <S(t)S(0)>
   	for(int step=0;step<datapoints;step++){
   			daten << step*delta << "\t" << sxsx_mean[step] << "\t";
   			daten << sysy_mean[step] << "\t";
   			daten << szsz_mean[step];
   			daten << "\t" << sx_mean[step] << "\t";
  			daten << sy_mean[step] << "\t";
			daten << sz_mean[step] << "\t";
			daten << sysx0[step];

			//Fehler
//			daten << "\t" << sqrt(sx_err[step]/(max_configs*(max_configs-1))) << "\t";
//			daten << sqrt(sy_err[step]/(max_configs*(max_configs-1))) << "\t";
//			daten << sqrt(sz_err[step]/(max_configs*(max_configs-1)));

  			daten << endl;
   	}

	delete[] y0;
	delete[] sx_mean;
	delete[] sxsx_mean;
	delete[] sxsx_err;
	delete[] sy_mean;
	delete[] sysy_mean;
	delete[] sysy_err;
	delete[] sz_mean;
	delete[] szsz_mean;
	delete[] szsz_err;
	delete[] sysx0;
}

int main(int argN, char** args){

	double tmax=100;
	int N=100;
	int ensembles=10e5;
	int B_function=1;
	double gamma=0.01;
	double H=1.;
	int datapoints=1000;
	string zusatz = "test";
	int coupling = 0;			//coupling =  	 1: equi	2: gamma	3: posneg
	double t_pulse=10.;
	if(argN == 11)
	{
		ensembles=atoi(args[1]);
		tmax=atoi(args[2]);
		N=atoi(args[3]);
		zusatz=args[4];
		coupling=atoi(args[5]);
		gamma=atof(args[6]);
		H=atof(args[7]); // H=0: Kein B-Feld, Feld in z-Richtung
		datapoints=atoi(args[8]); //Anzahl Datenpunkte bis t=tmax
		B_function=atoi(args[9]); //B-Feld function, 1: B=const
		t_pulse=atof(args[10]);
	}
	else{
		cout << "Falsche Anzahl Parameter. Verlangt: 10" << endl;
		return -1;
	}
	string name="RK_";
	name+=to_string(N);
	name+="spins_";
	name+=zusatz;
	name+=".txt";

	dynamic_RK(N,ensembles,tmax,name,gamma, coupling,H,datapoints,B_function, t_pulse);


	return 0;
}

