/*
 * RKtest_harmOszi.cpp
 *
 *  Created on: 13.04.2016
 *      Author: huedepohl
 */

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
#include <vector>
#include <chrono>
#include <random>
#include <iomanip>
#include <functional>
#include "../RK4_solver/dglsolver.h"
#include <cassert>
using namespace std;
//#include <iomanip>
#include <string>



void product(const vector<double>& a,const vector<double>& b, vector<double>&res){
	res[0]=a[1]*b[2]-a[2]*b[1];
	res[1]=a[2]*b[0]-a[0]*b[2];
	res[2]=a[0]*b[1]-a[1]*b[0];
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

void rightside(double w0,const vector<double>& x0, vector<double>& res){
	assert(x0.size()==2);
	assert(res.size()==2);
	res[0]=x0[1];
	res[1]=-w0*w0*x0[0];
}


int main(int argN, char** args){

	double tmax=10000;
	double h=0.1;
	double w=1.;
	int steps=tmax/h;
	vector<double> x;
	vector<double> res(steps,0.);
	vector<double> w0(2,0.),w1(2,0.),w2(2,0.),w3(2,0.);
	fstream daten("oszi_RKDavid.txt",ios::out);

	while(h>1e-3){
		x={0.,1.};
		for(unsigned i=0;i<steps;i++){
			if(i==steps-1) daten << h << "\t" << abs((x[0]-sin(w*h*i))/sin(w*h*i)) << endl;
			dglsolver::rk4step(h,x,x,w0,w1,w2,w3,bind(rightside,w,placeholders::_1,placeholders::_2));
		}
		h*=0.95;
	}

	return 0;
}






