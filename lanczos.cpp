/*
 * lanczos.cpp
 *
 *  Created on: 01.12.2015
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
#include <iomanip>
#include <string>
#include <eigen3/Eigen/Dense>
using namespace Eigen;
using Eigen::MatrixXd;

unsigned seed1 = (unsigned)chrono::system_clock::now().time_since_epoch().count();
//default_random_engine generator;
//random_device rd; //Seed fuer Mersenne-Twister
mt19937 gen(seed1); //Mersenne-Twister
//uniform_real_distribution<double> gleich(0, 1); //Gleichverteilung
normal_distribution<double> normal(0.0, 1./2.); //Gaussverteilung
//mu	sigma


vector<double> operator * (const double& faktor, const vector<double>& v1)
{
	vector<double> result = v1;
	for (unsigned int i = 0; i < result.size(); ++i)
	{
		result[i] *= faktor;
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


vector<double> operator + (const vector<double>& v1, const vector<double>& v2)
{
	vector<double> result = v1;
	for (unsigned int i = 0; i < result.size(); ++i)
	{
		result[i] += v2[i];
	}
	return result;
}

vector<double> operator - (const vector<double>& v1, const vector<double>& v2)
{
	vector<double> result = v1;
	for (unsigned int i = 0; i < result.size(); ++i)
	{
		result[i] -= v2[i];
	}
	return result;
}

double operator * (const vector<double>& v1, const vector<double>& v2)
{
	double result = 0;
	for (unsigned int i = 0; i < v1.size(); i++){
		result += v1[i] * v2[i];
	}
	return result;
}

vector<double> operator * (const MatrixXd &A, const vector<double>& v)
{	int size=v.size();
	vector<double> result(size,0.);
	for (int i = 0; i < size; i++){
		for (int j = 0; j < size; j++){
		result[i] += A(i,j) * v[j];
		}
	}
	return result;
}

double norm(vector<double> &v){
	int size=v.size();
	double result=0.;
	for(int i=0;i<size;i++){
		result+=v[i]*v[i];
	}
	return sqrt(result);
}

double fak(int k){
	if(k==0) return 1;
	if(k==1) return 1;
	else return k*fak(k-1);
}

MatrixXd Lanczos(const MatrixXd &A,int size,int steps){
	double numzero=1e-12; //lower limit for lanczos
	vector<vector<double>> v(steps+1,vector<double>(size,0.));
	vector<vector<double>> w(steps,vector<double>(size,0.));
	vector<vector<double>> w_(steps,vector<double>(size,0.));
	vector<double> alpha(steps,0.);
	vector<double> beta(steps,0.);
	MatrixXd tri(steps,steps);
	//random v1
	for(unsigned i=0;i<size;i++){
		v[1][i]=normal(gen);
	}
	v[1]=1./norm(v[1])*v[1];
	//lanczos algorithm
	cout << endl << "scalar products v_i*v_i+1 of Lanczos: " << endl;
	for(unsigned j=0;j<steps-1;j++){
		w_[j]=A*v[j+1];
		tri(j,j)=w_[j]*v[j+1];
		if(j==0) w[j]=w_[j]-tri(j,j)*v[j+1];
		else{ w[j]=w_[j]-tri(j,j)*v[j+1]-beta[j]*v[j];}
		beta[j+1]=norm(w[j]);
		if(fabs(beta[j+1])>numzero){
			v[j+2]=1./beta[j+1]*w[j];
			cout << v[j+1]*v[j+2] << endl;
		}
		else{
			break;
		}
	}
	cout << endl;
	w_[steps-1]=A*v[steps];
	tri(steps-1,steps-1)=w_[steps-1]*v[steps];
	//fill in rest of matrix
	for(unsigned i=0;i<steps;i++){
		for(unsigned j=0;j<steps;j++){
			if((i-j)*(i-j)==1) tri(i,j)=beta[max(i,j)];
		}
	}
	return tri;
}

vector<double> Poly_multi(const vector<double>& coeff1,const vector<double>& coeff2){
	int size1=coeff1.size();
	int size2=coeff2.size();
	vector<double> res(size1+size2-1,0.);
	for(unsigned i=0;i<size1;i++){
		for(unsigned j=0;j<size2;j++){
			res[i+j]+=coeff1[i]*coeff2[j];
		}
	}
	return res;
}

vector<double> Poly_mult_x(const vector<double> &coeff){
	int size=coeff.size();
	vector<double>res(size+1,0.);
	for(unsigned i=0;i<size;i++){
		res[i+1]=coeff[i];
	}
	return res;
}

//integrate polynom from -1 to 1
double integrate_Leg(const vector<double>& coeff){
	int size=coeff.size();
	double res=0;
	vector<double>int_coeff(size+1,0.);
	for(unsigned i=0;i<size;i++){
		int_coeff[i+1]=coeff[i]/(i+1);
//		res+=int_coeff[i+1];
//		res-=pow(-1,i+1)*int_coeff[i+1];
		if((i+1)%2!=0) res+=2*int_coeff[i+1];
	}
	return res;
}

//integrate polynom from -1/2 to 1/2
double integrate_GS(const vector<double>& coeff){
	int size=coeff.size();
	double res=0;
	vector<double>int_coeff(size+1,0.);
	for(unsigned i=0;i<size;i++){
		int_coeff[i+1]=coeff[i]/(i+1);
//		res+=int_coeff[i+1];
//		res-=pow(-1,i+1)*int_coeff[i+1];
		if((i+1)%2!=0) res+=2*(int_coeff[i+1]*pow(0.5,i+1));
	}
	return res;
}

//Koeffizienten der Legendre-Polynome als Vektor
vector<double> P_coeff_Leg(int n){
	vector<double> res(n+1,0.);
	double k=0;
	for(unsigned i=0;i<n+1;i++){
		if(n%2==0){
			if(i%2!=0){res[i]=0;}
			else{
				k=double(n-i)/2.;
				res[i]=pow(-1.,k)*fak(int(2*n-2*k))/(fak(n-k)*fak(n-2*k)*fak(k)*pow(2,n));
			}
		}
		if(n%2!=0){
			if(i%2==0){res[i]=0;}
			else{
				k=double(n-i)/2.;
				res[i]=pow(-1.,k)*fak(int(2*n-2*k))/(fak(n-k)*fak(n-2*k)*fak(k)*pow(2,n));
			}
		}
		res[i]*=sqrt((2.*n+1.)/2.);
	}
	return res;
}

//Koeffizienten der Polynome aus Gram-Schmidt (orthonormal)
vector<double> P_coeff_GS(int n){
	vector<double> res(n+1,0.);
	vector<double> tmp(n+1,0.);
	double int1=0;
	double int2=0;
	//Ausgangsbasis
	vector<vector<double>>w(n+1,vector<double>(n+1,0.));
	for(unsigned i=0;i<n+1;i++){
		w[i][i]=1.;
	}
	//Orthonormalbasis
	vector<vector<double>>v(n+1,vector<double>(n+1,0.));
	for(unsigned m=0;m<n+1;m++){
		tmp.assign(n+1,0.);
		for(unsigned i=0;i<m;i++){
			int1=integrate_GS(Poly_multi(v[i],w[m]));
			int2=integrate_GS(Poly_multi(v[i],v[i]));
			tmp=tmp+int1/int2*v[i];
		}
		v[m]=w[m]-tmp;
	}
	res=1./sqrt(integrate_GS(Poly_multi(v[n],v[n])))*v[n];
	return res;
}

int main(){
//	int size=1000; //size of matrix to tridiagonalize
//	int steps=10; // steps of lanczos
//	MatrixXd C(size,size);
//	MatrixXd Tri(size,size);
//	for(unsigned i=0;i<size;i++){
//		for(unsigned j=0;j<size;j++){
//			if(i==j) C(i,j)=i+1; //matrix to tridiagonalize
//		}
//	}
//	//SelfAdjointEigenSolver<MatrixXd> eigensolver(C);
//	//if(eigensolver.info()!= Success) abort();
//	Tri=Lanczos(C,size,steps);
//	cout << "Eigenwerte nach Eigen: " << endl;
//	for(int i=0;i<size;i++){
//		//cout << eigensolver.eigenvalues()[i] << endl;
//	}
//	cout << "Tridiagonalmatrix mit Lanczos:" << endl << Tri << endl;
//	SelfAdjointEigenSolver<MatrixXd> eigensolver2(Tri);
//	//if(eigensolver.info()!= Success) abort();
//	cout << "Eigenwerte der Tridiagonalmatrix:" << endl << eigensolver2.eigenvalues() << endl;

	int size=1;
		do{
		cout << endl << "Size of H: ";
		cin >> size;
		cout << endl;
		cout << "Calculate H: f --> x*f" << endl << "with orthonormal-Polynoms f" << endl;
		vector<double>P1(size,0.);
		vector<double>P2(size,0.);

		MatrixXd H(size,size);
		MatrixXd Tri(size,size);
		for(unsigned i=0;i<size;i++){
			for(unsigned j=0;j<size;j++){
				H(i,j)=integrate_GS(Poly_multi(P_coeff_GS(i),Poly_mult_x(P_coeff_GS(j))));
			}
		}
		Tri=Lanczos(H,size,size);

		SelfAdjointEigenSolver<MatrixXd> eigensolver(H);
		SelfAdjointEigenSolver<MatrixXd> eigensolver2(Tri);
		cout << "H:" << endl << H << endl << endl;
		cout << "Eigenvalues of H: " << endl;
		cout << eigensolver.eigenvalues() << endl;
		cout << endl;
		cout << "after Lanczos:" << endl << Tri << endl;
	}
	while(size!=0);


//	for(unsigned i=0;i<10;i++){
//		for(unsigned j=0;j<10;j++){
//			cout << endl;
////			if(i==j)
//			{
//				cout << i << "\t" << j << endl;
//				cout << integrate_GS(Poly_multi(P_coeff_GS(i),P_coeff_GS(j)));
//				cout << endl;
//			}
//		}
//	}

	return 0;
}

