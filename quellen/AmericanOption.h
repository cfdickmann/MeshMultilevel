#ifndef AMERICANOPTION_H_
#define AMERICANOPTION_H_

#include <stdio.h>
#include <vector>
#include <math.h>
#include <stdlib.h>
#include <time.h>
#include <iostream>
#include <fstream>
#include <cstring>

#include "RNG.h"
#include "Hilfsmittel.h"
#include "EuroBewerter.h"

#define MAX_CALL 1
#define MIN_PUT 0

namespace std {

class AmericanOption {
public:
	void run();
	void simplified();
	AmericanOption();
	virtual ~AmericanOption();
	vector<vector<double>> Levelergs;
	vector<vector<double>> LevelergsFiner;

	int option; // MAX_CALL or MIN_PUT
	double delta; //dividend yield
	vector<double> X0; // Spot
	double Strike; // Ausuebungspreis
	double r; // interest rate
	vector<double> sigma; //Volatility
	double T; //Gesamtzeit
	int N; //time discretization
	int D;
	vector<vector<vector<vector<double>>>> X;
	vector<vector<vector<double>>> V;
	void Daten();

	double dt;
	int L;
	vector<int> n;

	double Pfad(vector<vector<double>> X, int l);
	bool Kernel(vector<double> x, vector<double> y, int D, double threshold);
	void Pfadgenerieren(vector<vector<double>> X, int start, vector<double> S, RNG* generator);
	void Pfadgenerieren(vector<vector<double>> X, vector<vector<double>> wdiff, int start, vector<double> S);
	void printInfo();
	void addLevelPath(int l);
	void ErgebnisseAusgeben(int l);
	vector<vector<vector<double>>> weights;
	vector<vector<vector<double>>>  weight_sum;
	void weights_erstellen(int l);
	EuroBewerter EB;
	int Mtraining(int l);
	double kernelD(vector<double>von, vector<double> nach, double dt);
	double payoff(vector<vector<double>> x, int time);
	double payoff(vector<double> x, int time);
	void trainingpaths_erstellen(int l);
	double C_estimate_Mesh(vector<double> x, int lauf, int l);
	void trainingpaths_regression(int l);
};

}
#endif
