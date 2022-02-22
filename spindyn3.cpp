/*
 * spindyn3.cpp
 *
 *  Created on: 16.10.2015
 *      Author: huedepohl
 */

//Berechnet die Dynamik von Zentralspin und aller Badspins einzeln durch Auswertung der
//analytischen Funktion S(t)

#include <iostream>
#include <fstream>
#define _USE_MATH_DEFINES
#include <math.h>
#include <vector>
#include <chrono>
#include <random>
using namespace std;
#include <string>

unsigned seed1 = chrono::system_clock::now().time_since_epoch().count();
	//default_random_engine generator;
	random_device rd; //Seed fuer Mersenne-Twister
	mt19937 gen(seed1); //Mersenne-Twister
	//uniform_real_distribution<double> gleich(0, 1); //Gleichverteilung
	normal_distribution<double> normal(0.0, 1./2.); //Gaussverteilung
									  //mu	sigma

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

vector<double> operator + (const vector<double>& v1, const vector<double>& v2)
{
	vector<double> result=v1;
	for (unsigned int i = 0; i < result.size(); ++i)
	{
		result[i] += v2[i];
	}
	return result;
}

double operator * (const vector<double>& v1, const vector<double>& v2)
{
	double result=0;
	for(unsigned int i=0;i<3;i++){
		result+=v1[i]*v2[i];
	}
	return result;
}

void out(vector<double>& v){
	cout << v[0] << endl << v[1] << endl << v[2] << endl;
}

void S(double t,double omega,const vector<double>& S0, const vector<double>& n, vector<double>& res){
	res=(S0*n)*n;
	res=res+cos(omega*t)*(S0+(-1)*(S0*n)*n);
	res=res+sin(omega*t)*product((S0+(-1)*(S0*n)*n),n);
}

void dynamic_direct(int N,int max_configs, double h, int max_step, string name){
	fstream daten(name,ios::out);
	daten << "#t \t <Sz(t)Sz(0)>" << endl;
	//save CS dyn for mean
	vector<double> Sz_mean(max_step+1,0.);
	//initial condition
	//Bath-Spins
	vector<vector<double>> S(N+1,vector<double>(3,0.));
	vector<vector<double>> S0(N+1,vector<double>(3,0.));
	double Jq=1.;
	//couplings
	vector<double> J(N,0.);
	//overhauser
	vector<double> A(3,0.);
	double mu_B=1.;
	double g_e=1.;
	double omega_cs=0.;
	double omega_bath=0.;
	vector<double> n_A(3,0.);
	vector<double> n_CS(3,0.);
	string ausgabe="Fortschritt ";
	ausgabe+=name;
	ausgabe+= "bei: ";

	//Initialize couplings
	for(int i=0;i<N;i++){
		J[i]=sqrt(6.*N/(2.*N*N+3.*N+1.))*(N-i)/N*Jq;
	}

	//Calculate Dynamic
	for(int config=0;config<max_configs;config++){
		if(config%(max_configs/10)==0){
			cout << ausgabe << int(config*10/max_configs) << "/10" << endl;
		}
		//Initialize random Config for all spins
		for(int n=0;n<N+1;n++){
			for(int comp=0;comp<3;comp++){
				S[n][comp]=normal(gen);
			}
		}
		//starting values
		S0=S;
		for(int step=0;step<max_step;step++){
			Sz_mean[step]+=S0[0][2]*S[0][2]/max_configs;
			A.assign(3,0.);
			//Calculate Overhauser-Field
			for(int i=1;i<N+1;i++){
				A=A+J[i]*S[i+1];
			}
			//unit vector in direction of overhauser
			n_A=1./sqrt(A*A)*A;
			//unit vector in direction of CS
			n_CS=1./sqrt(S[0]*S[0])*S[0];
			//Larmor frequency of CS
			omega_cs=mu_B*g_e*sqrt(A*A);
			//Larmor frequency of BSs
			omega_bath=mu_B*g_e*sqrt(S[0]*S[0]);
			//update CS
			S(step*h,omega_cs,S0[0],n_A,S[0]);
			//update bath spins
			for(int i=1;i<N+1;i++){
				S(step*h,omega_bath,S0[i],n_CS,S[i]);
			}
		}
	}
	for(int step=0;step<max_step;step++){
		daten << step*h << "\t" << Sz_mean[step] << endl;
	}
}

int main(int argN, char** args){
	int N = 10;
	int ensembles=10e5;
	int steps=1000;
	string zusatz = "test";

	if(argN == 5)
	{
		N = atoi(args[1]);
		ensembles=atoi(args[2]);
		steps=atoi(args[3]);
		zusatz=args[4];

	}
	else{
		cout << "Zu wenig Parameter. Verlangt: 4" << endl;
		return -1;
	}

	string dateiname = "";
	dateiname += to_string(N);
	dateiname += "spins";
	dateiname += zusatz;
	dateiname += ".txt";
	cout << endl << "Rechnung gestartet mit" << endl;
	cout << "Spins: " << N << endl;
	cout << "Mittelung über " << ensembles << " Werte" << endl;
	cout << "Zeitschritte: " << steps << endl;
	cout << "Dateiname: " << dateiname << endl;
	dynamic(N, ensembles, 0.01,steps ,dateiname);

	return 0;
}








