/*
 * start_value_gen.cpp
 *
 *  Created on: 12.04.2016
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
using namespace std;
//#include <iomanip>
#include <string>

unsigned seed1 = chrono::system_clock::now().time_since_epoch().count();
	//default_random_engine generator;
	//random_device rd; //Seed fuer Mersenne-Twister
	mt19937 gen(seed1); //Mersenne-Twister
	//uniform_real_distribution<double> gleich(0, 1); //Gleichverteilung
	normal_distribution<double> normal(0.0, 1./2.); //Gaussverteilung



int main(){
	int N;
	double zahl=0.;
	cout << "Anzahl an Gaußverteilten Zufallszahlen mit Mittelwert 0 und Varianz 1/4:" << endl;
	cin >> N;
	fstream daten("rand.txt",ios::out);
	for(int i=0;i<N;i++){
		daten << normal(gen) << endl;
	}

	return 0;
}
