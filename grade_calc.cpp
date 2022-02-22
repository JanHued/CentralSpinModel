/*
 * gradecalc.cpp
 *
 *  Created on: 26.01.2016
 *      Author: huedepohl
 */

#include <float.h>

#include <iostream>
#include <fstream>
#define _USE_MATH_DEFINES
#include <math.h>
using namespace std;
#include <iomanip>
#include <string>
#include <vector>


void ask_for_values(){
	cout << "Insert Credits of Module and Grade as:" << endl << "credits" << endl << "grade" << endl;
}

bool legit_cred(double credits){
	bool val=credits>0;
	val*=credits<=60;
	return val;
}

bool legit_grade(double grade){
	bool val = grade >=1.0;
	val*=grade<=4.0;
	return val;
}

void tell_score(int score){
	cout << "Current Credit Score: " << score << endl;
}

void tell_grade(double grade){
	cout << "Final Grade:" << endl;
	cout << "-------- " << grade << " --------" << endl;
}

int main(){
	cout << "%%%%% GradeCalc 2016 %%%%%" << endl << endl;
	vector<vector<double> > modules;
	int score=0;
	int max_credits=120;
	double av_grade=0;
	vector<double> temp(2,0.);		// Credits	Grade
	ask_for_values();
	while(score<max_credits){
		while(not(legit_grade(temp[1])&&legit_cred(temp[0]))){
			cin >> temp[0];
			cin >> temp[1];
		}
		modules.push_back(temp);
		score+=temp[0];
		temp.assign(2,0.);
		tell_score(score);
	}
	for(int i=0;i<modules.size();i++){
		av_grade+=modules[i][0]*modules[i][1];
	}
	av_grade/=score;
	tell_grade(av_grade);
	return 0;
}


