#ifndef AMERICANOPTION_H_
#define AMERICANOPTION_H_

#include <iostream>
#include <stdio.h>
#include <stdio.h>
#include <vector>
#include <stdio.h>
#include <math.h>
#include <iostream>
#include <stdlib.h>
#include <time.h>
#include "math.h"
#include "Hilfsmittel.h"
#include <stdlib.h>
#include <stdio.h>
#include <math.h>
#include <math.h>
#include <iostream>
#include <fstream>
#include <stdlib.h>
#include <cstring>
#include "RNG.h"

#include "EuroBewerter.h"

#define MAX_CALL 1
#define MIN_PUT 0
#define ITO 1
#define ITOrho 12
#define EULER 2
#define LIBOR 5
namespace std {

class AmericanOption {
public:
	void run();
	AmericanOption();
	virtual ~AmericanOption();
	vector<double>* Levelergs;
//	vector<double>* cv;

	void printInfo();
	void addLevelPath(int l);

	void ErgebnisseAusgeben(int l);
	double getLevelVar(int l);
	double getLevelComp(int l);

	int option; // MAX_CALL or MIN_PUT
	double delta; //dividend yield
	double* X0; // Spot
	double Strike; // Ausuebungspreis
	double r; // interest rate
	double* sigma; //Volatility
	double T; //Gesamtzeit

	double**** X;
	int N; //time discretization
	int D;
	double dt;
	int L;
int * n;

	double Pfad(double ** X, int l);
	bool Kernel(double *x, double* y, int D, double threshold);
	void Pfadgenerieren(double** X, int start, double* S, RNG* generator);
	void Pfadgenerieren(double** X, double** wdiff, int start, double * S);
	void Daten();
	double *** weights;
	double *** weight_sum;
	void weights_erstellen(int l);
	EuroBewerter EB;

	double max(double d1, double d2);

	void simplified();
	double payoff(double** x, int time);
	double payoff(double* x, int time);
	void trainingpaths_erstellen(int l);
	double C_estimate_Mesh(double* x, int lauf, int l);
	void trainingpaths_regression(int l);
	double ***V;
	int Mtraining(int l);
	double kernelD(double * von, double* nach, double dt);
};

} /* namespace std */
#endif /* AMERICANOPTION_H_ */
