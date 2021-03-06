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
	vector<double>* Levelergs;
	vector<double>* LevelergsFiner;

	int option; // MAX_CALL or MIN_PUT
	double delta; //dividend yield
	double* X0; // Spot
	double Strike; // Ausuebungspreis
	double r; // interest rate
	double* sigma; //Volatility
	double T; //Gesamtzeit
	int N; //time discretization
	int D;
	double**** X;
	double ***V;
	void Daten();

	double dt;
	int L;
	int * n;

	double Pfad(double ** X, int l);
	bool Kernel(double *x, double* y, int D, double threshold);
	void Pfadgenerieren(double** X, int start, double* S, RNG* generator);
	void Pfadgenerieren(double** X, double** wdiff, int start, double * S);
	void printInfo();
	void addLevelPath(int l);
	void ErgebnisseAusgeben(int l);
	double *** weights;
	double *** weight_sum;
	void weights_erstellen(int l);
	EuroBewerter EB;
	int Mtraining(int l);
	double kernelD(double * von, double* nach, double dt);
	double payoff(double** x, int time);
	double payoff(double* x, int time);
	void trainingpaths_erstellen(int l);
	double C_estimate_Mesh(double* x, int lauf, int l);
	void trainingpaths_regression(int l);
};

}
#endif
