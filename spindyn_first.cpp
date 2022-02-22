/*
 * spindyn_first.cpp
 *
 *  Created on: 01.10.2015
 *      Author: huedepohl
 */
#include <iostream>
#include <fstream>
#define _USE_MATH_DEFINES
#include <math.h>
#include <vector>
#include <chrono>
#include <random>
using namespace std;
//#include <iomanip>
#include <string>

unsigned seed1 = chrono::system_clock::now().time_since_epoch().count();
	//default_random_engine generator;
	//random_device rd; //Seed fuer Mersenne-Twister
	mt19937 gen(seed1); //Mersenne-Twister
	//uniform_real_distribution<double> gleich(0, 1); //Gleichverteilung aus Mersenne
	normal_distribution<double> normal(0.0, 1.0); //Gaussverteilung aus Mersenne

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

void out(vector<double> v){
	cout << v[0] << endl << v[1] << endl << v[2] << endl;
}

void dynamic(int N, double h, int max_step, string name){
	fstream daten(name,ios::out);
	daten << "#t \t x \t y \t z \t E" << endl;
	vector<vector<double> > s(max_step,vector<double>(3,0.));
	vector<vector<double> > s_mean(max_step,vector<double>(3,0.));
	vector<vector<double> > s0(max_step,vector<double>(3,0.));
	vector<double> energy(max_step,0.);
	vector<double> s_start={0,0,1}; //Up_init
	vector<double> B(3,0.);
	//vector<double> sum(3,0.);
	vector<double> k1(3,0.);
	vector<double> k2(3,0.);
	vector<double> k3(3,0.);
	vector<double> k4(3,0.);
	double mu=0;
	double sigma2=0;

	//Initial Condition
	s0[0]=s_start;
	//Calculate Dynamic
	for(int i=0;i<N;i++){
		if(i%(N/10)==0){
			cout << "Fortschritt: " << int(i*10/N) << "/10" << endl;
		}
		//Initialize
		s=s0;
		for(int comp=0;comp<3;comp++){
			B[comp]=normal(gen);
			mu+=B[comp];
			sigma2+=B[comp]*B[comp];
		}
		for(int step=0;step<max_step-1;step++){
			//RK-4
			k1=h*product(B,s[step]);
			k2=h*product(B,s[step]+1/2*k1);
			k3=h*product(B,s[step]+1/2*k2);
			k4=h*product(B,s[step]+k3);
			s[step+1]=s[step]+1/6.*(k1+2*k2+2*k3+k4);
			s_mean[step]=s_mean[step]+1./N*s[step];
			energy[step]+=1./N*s[step]*B;
		}
	}
	mu/=3*N;
	sigma2/=3*N;
	for(int step=0;step<max_step-1;step++){
		daten << h*step << "\t" << s_mean[step][0] << "\t" << s_mean[step][1] << "\t" << s_mean[step][2] << "\t" << energy[step] << endl;
	}
	cout << "mu: " << mu << endl;
	cout << "sigma^2: " << sigma2 << endl;
}

int main(){
	dynamic(1000,0.01,1000,"spindyn_s0_1e3.txt");
	dynamic(10000,0.01,1000,"spindyn_s0_1e4.txt");
	dynamic(100000,0.01,1000,"spindyn_s0_1e5.txt");



	return 0;
}



