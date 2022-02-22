/*
 * hierarchy_density.cpp
 *
 *  Created on: 17.02.2016
 *      Author: huedepohl
 */

#include <float.h>
#define DEBUG 0
#include <iostream>
#include <fstream>
#define _USE_MATH_DEFINES
#include <math.h>
#include <vector>
#include <chrono>
#include <random>
using namespace std;
#include <iomanip>
#include <string>
#include <eigen3/Eigen/Dense>
#include <functional>
#include <cassert>
#include <gsl/gsl_odeiv2.h>
#include <gsl/gsl_errno.h>

unsigned seed1 = (unsigned)chrono::system_clock::now().time_since_epoch().count();
//default_random_engine generator;
//random_device rd; //Seed fuer Mersenne-Twister
mt19937 gen(seed1); //Mersenne-Twister
//uniform_real_distribution<double> gleich(0, 1); //Gleichverteilung
normal_distribution<double> normal(0.0, 1./2.); //Gaussverteilung
//mu	sigma

//parameter of right side of ODEs
struct parameter{
	double* epsilon;
	int size;
	double* weight;
	double B;
	int B_function;
};

void out(const double* v, int v_size){
	for(int i=0;i<v_size;i++){
		cout << v[i] << endl;
	}
}

double Bfield(const double t,const double B0, const int function){
	assert(function>0);
	switch(function){
		case 1: return  B0; //B=const
	}
}

void couplings_gamma(double* J, const int J_size, const double gamma){
	assert(gamma<1);
	const int N=J_size;
	assert(N>1);
	const double A0=sqrt((1.-exp(-2.*gamma))/(1.-exp(-2.*gamma*N)));
	for(int i=0;i<N;i++){
		J[i]=A0*exp(-double(i)*gamma);
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

void couplings_equidist(double* J, const int J_size){
	int N=J_size;
	double Jq=1.;
	for(int i=0;i<N;i++){
			J[i]=sqrt(6.*N/(2.*N*N+3.*N+1.))*(N-i)/N*Jq;
	}
}

//evaluate polynom p at x with horner scheme
double poly_val(const double* p, const int p_size, const double x){
	if(p_size==2) return p[0]+p[1]*x;
	double res=x*p[p_size-1];
	for(int i=p_size-2;i>=1;i--){
		res+=p[i];
		res*=x;
	}
	res+=p[0];
	return res;
}

//vectorproduct for vectors of arbitrary length with given indices
inline void crossprod(const double* v1, const double coeff, const double* v2, double* dest,const int ind1, const int ind2, const int inddest){
	dest[inddest]=coeff*(v1[ind1+1]*v2[ind2+2]-v1[ind1+2]*v2[ind2+1]);
	dest[inddest+1]=coeff*(v1[ind1+2]*v2[ind2]-v1[ind1]*v2[ind2+2]);
	dest[inddest+2]=coeff*(v1[ind1]*v2[ind2+1]-v1[ind1+1]*v2[ind2]);
}

inline void crossprod(const double* v1, const double* v2, double* dest,const int ind1, const int ind2, const int inddest){
	dest[inddest]=v1[ind1+1]*v2[ind2+2]-v1[ind1+2]*v2[ind2+1];
	dest[inddest+1]=v1[ind1+2]*v2[ind2]-v1[ind1]*v2[ind2+2];
	dest[inddest+2]=v1[ind1]*v2[ind2+1]-v1[ind1+1]*v2[ind2];
}

void integrate_polynom(const double* polycoeff, int polycoeff_size, double* res){
	res[0]=0.;
	for(unsigned i=0;i<polycoeff_size;++i){
			res[i+1]=polycoeff[i]/(i+1);
	}
}

void differentiate_polynom(const double* polycoeff, const int polycoeff_size, double* res){
	for(unsigned i=0;i<polycoeff_size-1;++i){
			res[i]=polycoeff[i+1]*(i+1);
	}
}

double Newton(const double* poly, const int poly_size, const double* diff, const int diff_size, const double start){
	double err=1e-10;
	double x,xnext=start;
	while(abs(poly_val(poly, poly_size, xnext))>err){
		x=xnext;
		xnext=x-poly_val(poly,poly_size,x)/poly_val(diff,diff_size,x);
	}
	return xnext;
}

double calc_median(double* intcoeff, const int intcoeff_size, const double* polycoeff, const int polycoeff_size, const double a,const double b){
	double res=0;
	double x1=a;
	double x2=b;
	double* newton_poly=intcoeff;
	res=poly_val(intcoeff,intcoeff_size,b)+poly_val(intcoeff,intcoeff_size,a);
	res/=2;
	newton_poly[0]-=res;
	return Newton(newton_poly,intcoeff_size,polycoeff, polycoeff_size, (a+b)/2.);
}

double calc_llim(const double ulim, const double weight, const double* coeff, const int coeff_size){
	double llim=0.;
	double* intcoeff=new double[coeff_size+1];
	integrate_polynom(coeff,coeff_size, intcoeff);
	double F_x=poly_val(intcoeff,coeff_size+1, ulim);
	llim=1./(intcoeff[coeff_size])*(F_x-weight);
	llim=pow(llim,1./(coeff_size));
	return llim;
}

//calculates parameters for continous spectral density
void calc_param_continous(const int coupling, const int N, const double gamma, const int tmax, const bool expdisc, const int intervals, double faktor, double* epsilon, double* w){

	double* limits=new double[intervals+1];
	double interval=0.;
	double* dens;
	int dens_size;
	double a_max;
	double a_min=0.;
	//boundaries of integral
	double llim,ulim=0.;

	//set min and max of coupling and related density-polynom
	switch(coupling){
		case 1: a_max=sqrt(3./N);
				dens_size=3;
				dens=new double[dens_size];
				dens[0]=0.;
				dens[1]=0.;
				dens[2]=pow(N,3./2.)/sqrt(3);
				break;
		case 2: a_max=sqrt(2*gamma); //N=infty
				dens_size=2;
				dens=new double[dens_size];
				dens[0]=0.;
				dens[1]=1./gamma;
				break;
	}

	double* int_dens=new double[dens_size+1];
	integrate_polynom(dens,dens_size,int_dens);
	ulim=a_max;
	interval=a_max-a_min;
	double intwidth=interval/intervals;
	//Zunächst äquidistante Diskretisierung
	for(unsigned i=0;i<intervals;++i){
		limits[i]=a_max-i*intwidth;
	}
	limits[intervals]=0.;
	if(expdisc){
			if(faktor<1){
				for(unsigned i=0;i<=intervals;i++){
					limits[i]*=pow(faktor,i);
				}
			}
			else if(faktor==1){
				//faktor an lvl anpassen, sodass Rechnung "korrekt" bis tmax
				faktor=pow(1./tmax/intwidth,1./intervals);
				for(unsigned i=0;i<=intervals;i++){
					limits[i]*=pow(faktor,i);
				}
			}
		}
	//Koeffizienten berechnen
	for(unsigned i=0;i<intervals;i++){
		w[i]=poly_val(int_dens,dens_size+1,limits[i])-poly_val(int_dens,dens_size+1, limits[i+1]);
		epsilon[i]=calc_median(int_dens,dens_size+1,dens,dens_size,limits[i+1],limits[i]);
	}
	delete[] limits;
	delete[] dens;
	delete[] int_dens;
}

//calculates parameters for discrete spectral density
void calc_param_discrete(const int coupling, const int N, const double gamma, const int tmax, const bool expdisc, const int intervals, double faktor, double* epsilon, double* w){
	double* J=new double[N];
	switch(coupling){
		case 1: couplings_equidist(J,N);break;
		case 2: couplings_gamma(J,N,gamma);break;
	}
	double* limits=new double[intervals+1];
	double interval=J[0];
	double intwidth=interval/intervals;
	for(unsigned i=0;i<=intervals;++i){
		limits[i]=interval-i*intwidth;
	}
	if(expdisc){
		if(faktor<1){
			for(unsigned i=0;i<=intervals;i++){
				limits[i]*=pow(faktor,i);
			}
		}
		else if(faktor==1){
			//faktor an lvl anpassen, sodass Rechnung gut bis tmax
			faktor=pow(1./tmax/intwidth,1./intervals);
			for(unsigned i=0;i<=intervals;i++){
				limits[i]*=pow(faktor,i);
			}
		}
	}
	int index=0;
	double J2=0.;
	double J3_sum=0.;
	for(unsigned i=0;i<intervals;i++){
		while(J[index]>limits[i+1]){
			J2=J[index]*J[index];
			w[i]+=J2;
			J3_sum+=J2*J[index];
			index++;
		}
		epsilon[i]=J3_sum/w[i];
		J3_sum=0.;
	}
	delete[] J;
	delete[] limits;
}

int rightside(double t, const double* y_0, double* dest, void* params){
	parameter* p=(parameter*) params;
	const double* epsilon=p->epsilon;
	const int lvl=p->size;
	const double* weight=p->weight;
	double H=p->B;
	const int B_function=p->B_function;
	H=Bfield(t,H,B_function);

	switch(lvl){
		//Nur für 0 Intervalle ("Merkulov-Fall")
		case 0:	crossprod(y_0,y_0,dest,3,0,0);
				dest[3]=0;
				dest[4]=0;
				dest[5]=0;
				break;
		default:
			double* P1=new double[3];
			P1[0]=0.;
			P1[1]=0.;
			P1[2]=0.;
			for(int i=0;i<lvl;i++){
				P1[0]+=sqrt(weight[i])*y_0[3*(i+1)];
				P1[1]+=sqrt(weight[i])*y_0[3*(i+1)+1];
				P1[2]+=sqrt(weight[i])*y_0[3*(i+1)+2];
			}
			//Dynamik von S0 berechnen
			crossprod(P1,y_0,dest,0,0,0);
			//Dynamik durch B-Feld
			dest[0]+=H*y_0[1];
			dest[1]-=H*y_0[0];
			//Dynamik der Q_i berechnen
			for(int i=1;i<=lvl;++i){
				crossprod(y_0,epsilon[i-1],y_0,dest,0,3*i,3*i);
			}
			delete[] P1;
	}

	return GSL_SUCCESS;
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

void solveODEs(const int N,const int intervals,double faktor, const bool contdens, const int max_configs, const int tmax, const int coupling, const float gamma, const string name, const int datapoints, const double H, const int B_function, const double t_pulse){

	const double delta=double(tmax)/datapoints;
	double t,ti;

	fstream daten(name, ios::out);

	double* epsilon=new double[intervals];
	double* w=new double[intervals];
	bool expdisc=true;

	if(contdens){
		calc_param_continous(coupling, N, gamma, tmax, expdisc, intervals, faktor, epsilon, w);
	}
	else{
		calc_param_discrete(coupling, N, gamma, tmax, expdisc, intervals, faktor, epsilon, w);
	}


	double weight=0.;
	cout << "Median epsilon: " << "\t" << "Gewicht W:" << endl;
	for(int i=0;i<intervals;i++){
		weight+=w[i];
		cout << epsilon[i] << "\t\t" << w[i] << endl;
	}
	cout << "Gesamtgewicht: " << weight << endl << endl;

	//P-vectors with Q[0,1,2]=CS
	unsigned dim=intervals ? (3*(intervals+1)) : 6;
	double* y0=new double[dim];
	//starting value of S0_z
	double* S0=new double[3];

	//steps befor pulse
	int pulsestep=int(t_pulse/delta);

	// ensemble mean
	double* sx_mean=new double[datapoints];
	double* sxsx_mean=new double[datapoints];
	double sx_oldmean;
	double* sx_err=new double[datapoints];
	double* sy_mean=new double[datapoints];
	double* sysy_mean=new double[datapoints];
	double sy_oldmean;
	double* sy_err=new double[datapoints];
	double *sz_mean=new double[datapoints];
	double* szsz_mean=new double[datapoints];
	double sz_oldmean;
	double* sz_err=new double[datapoints];
	double* sysx0=new double[datapoints];


	double sxsx=0.;
	double sysy=0.;
	double szsz=0.;

	for(unsigned i=0;i<datapoints;++i){
		sx_mean[i]=0.;
		sxsx_mean[i]=0.;
		sx_err[i]=0.;
		sy_mean[i]=0.;
		sysy_mean[i]=0.;
		sy_err[i]=0.;
		sz_mean[i]=0.;
		szsz_mean[i]=0.;
		sz_err[i]=0.;
		sysx0[i]=0.;
	}

	parameter params;
	params.epsilon=epsilon;
	params.weight=w;
	params.size=intervals;
	params.B=H;
	params.B_function=B_function;

	parameter* p = &params;

	//initialize odes
	gsl_odeiv2_system ode = {rightside, NULL, dim, p};

	//rk4 = classic 4th order RK, rkf45 = explicit embedded Runge-Kutta-Fehlberg (4, 5) method.
	gsl_odeiv2_driver* d = gsl_odeiv2_driver_alloc_y_new(&ode,gsl_odeiv2_step_rkf45,1e-6,1e-6,0);

//	gsl_odeiv2_driver_set_hmin(d,1e-6);
//	gsl_odeiv2_driver_set_hmax(d,delta);
//	gsl_odeiv2_driver* d = gsl_odeiv2_driver_alloc_y_new(&ode,gsl_odeiv2_step_rk4,0.002,1e-6,0.)


	//ensemble average
	for (int config = 1; config <= max_configs; config++){
		t=0.;
		ti=0.;

		if(config%(max_configs/100)==0){
			cout.flush();
			cout << "\r" << "Fortschritt: " << config/(max_configs/100) << "%";
		}

//starting values for P-fields

		//random
		for(unsigned i=0;i<dim;i++){
				y0[i]=normal(gen);
		}


		flip_spin(y0,0);


		//starting values for correlation
		S0[0]=y0[0];
		S0[1]=y0[1];
		S0[2]=y0[2];

		for(int step=0;step<datapoints;step++){
			ti=(step+1)*delta;

			if(pulsestep){
				if(step%pulsestep==0){
					//orient S0 in x-dir
					flip_spin(y0,0);
				}
			}

			sxsx=y0[0]*S0[0];
			sysy=y0[1]*S0[1];
			szsz=y0[2]*S0[2];

			//S_0 mitteln ohne Fehler
			sx_mean[step]+=y0[0]/max_configs;
			sy_mean[step]+=y0[1]/max_configs;
			sz_mean[step]+=y0[2]/max_configs;

			sxsx_mean[step]+=sxsx/max_configs;
			sysy_mean[step]+=sysy/max_configs;
			szsz_mean[step]+=szsz/max_configs;

			sysx0[step]+=y0[1]*S0[0]/max_configs;

			//S_0 mitteln mit Fehler
//			sx_oldmean=sxsx_mean[step];
//			sxsx_mean[step]+=(sxsx-sxsx_mean[step])/config;
//			sx_err[step]+=(sxsx-sxsx_mean[step])*(sxsx-sx_oldmean);
//			sy_oldmean=sysy_mean[step];
//			sysy_mean[step]+=(sysy-sysy_mean[step])/config;
//			sy_err[step]+=(sysy-sysy_mean[step])*(sysy-sy_oldmean);
//			sz_oldmean=szsz_mean[step];
//			szsz_mean[step]+=(szsz-szsz_mean[step])/config;
//			sz_err[step]+=(szsz-szsz_mean[step])*(szsz-sz_oldmean);


            //time evolution
            int status = gsl_odeiv2_driver_apply(d,&t,ti,y0);

            if(status != GSL_SUCCESS) {
                cerr << "GSL Error: return value = " << status << endl;
                exit(1);
            }
		}
	}

	// 	save data of <Sz(t)Sz(0)>
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

    gsl_odeiv2_driver_free(d);

	delete[] y0;
	delete[] S0;
	delete[] epsilon;
	delete[] w;
	delete[] sx_mean;
	delete[] sx_err;
	delete[] sy_mean;
	delete[] sy_err;
	delete[] sz_mean;
	delete[] sz_err;
}

int main(int argN, char** args){

	int N=100;
	int ensembles=10e5;
	double gamma=0.01;
	double H=1.;
	int lvl=3;
	int tmax=1000;
	int coupling=2;
	bool denstype=0;
	double faktor=0.9;
	int datapoints=1000;
	int B_function=1;
	double t_pulse=10.;
	string zusatz = "test";

	if(argN == 14)
	{
		ensembles=atoi(args[1]);
		lvl=atoi(args[2]);
		tmax=atoi(args[3]);
		gamma=atof(args[4]);
		N=atoi(args[5]); 			// 1 := infty
		coupling=atoi(args[6]);     // 1:=equi	2:=gamma
		denstype=atoi(args[7]);		// 0: discrete	1: continous
		faktor=atof(args[8]);		// faktor intervallgrenzen spektraldichte 1: an lvl angepasst (sonst 0<faktor<1)
		zusatz=args[9];
		datapoints=atoi(args[10]);	//datapoints between t=0 and t=tmax
		H=atof(args[11]);				//amplitude of B
		B_function=atoi(args[12]); //B(t)	1: B=const
		t_pulse=atof(args[13]);		//0: no pulses
	}
	else{
		cout << "Falsche Anzahl Parameter. Verlangt: 13" << endl;
		return -1;
	}

	string dateiname="density_";
	if(coupling==1) dateiname+= "equi";
	else if(coupling==2) dateiname+= "gamma=" + string(args[4]);
	dateiname += "_intervals=";
	dateiname += to_string(lvl);
	dateiname += "_";
	if(N==1) dateiname+="INFTY";
	else dateiname += to_string(N);
	dateiname += "spins_";
	if(denstype==0) dateiname += "discdens_";
	else if(denstype==1) dateiname += "contdens_";
	dateiname+=zusatz;
	dateiname += ".txt";

	cout << "Calculation with density started with: " << endl;
	cout << ensembles << " Ensembles" << endl;
	cout << "Gamma=" << gamma << endl;
	cout << "intervals=" << lvl << endl;
	cout << "Daten schreiben in: " << dateiname << endl;

	solveODEs(N,lvl, faktor, denstype,ensembles, tmax,coupling, gamma, dateiname,datapoints,H,B_function,t_pulse);
	cout << endl << "---Finished---" << endl;

	return 0;
}
