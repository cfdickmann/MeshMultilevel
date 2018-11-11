#include "EuroBewerter.h"
#include "math.h"
#include <stdio.h>
#include <stdlib.h>
#define pi 3.141592654

namespace std {

double cnd(double x) {
	return 0.5 * (1 + erf(x / sqrt(2)));
}

double dnd(double x) {
	return 1. / sqrt(2. * 3.14159265) * exp(-x * x / 2.);
}

double EuroBewerter::call(double t, double T, double X0, double Strike,
		double r, double delta, double sigma) {
	double rStrich = r - delta;
	double TStrich = T - t;
	double d1 = (log(X0 / Strike) + (rStrich + sigma * sigma / 2.) * TStrich) / (sigma * sqrt(TStrich));
	double d2 = d1 - sigma * sqrt(TStrich);
	double wert = (X0 * cnd(d1) - Strike * exp(-rStrich * TStrich) * cnd(d2));
	return wert *exp(-(delta) * T)* exp(-t * (r - delta));
}

double EuroBewerter::put(double t, double T, double X0, double Strike, double r,
		double delta, double sigma) {
	double rStrich = r - delta;
	double TStrich = T - t;
	double d1 = (log(X0 / Strike) + (rStrich + sigma * sigma / 2) * TStrich)/ (sigma * sqrt(TStrich));
	double d2 = d1 - sigma * sqrt(TStrich);
	return exp(-r * t)	* (exp(-delta * TStrich)* (Strike * exp(-rStrich * TStrich) * cnd(-d2)	- X0 * cnd(-d1)));
}

//Margrabes Formel
double EuroBewerter::exchange_option(double x, double y, double t, double T,
		double r, double delta, double sigma) {
	// y to  be replaced
	double TStrich = T - t;

	double d1 = (log(x / y) + (sigma * sigma) * TStrich)
			/ (sqrt(2.) * sigma * sqrt(TStrich));
	double d2 = (log(x / y) - (sigma * sigma) * TStrich)
			/ (sqrt(2.) * sigma * sqrt(TStrich));

	double erg = x * exp(-delta * TStrich) * cnd(d1)
			- y * exp(-delta * TStrich) * cnd(d2);
	return erg * exp(-r * t);
}

double EuroBewerter::european_MaxCall_ND(double* x, int D, double t, double T,
		double Strike, double r, double delta, double sigma, double dt) {
	
	if(D==1)
		return call(t,T,x[0],Strike,r,delta,sigma);
	
	if(dt==0)
		return 0;
		
	double d_minus[D];
	double d_plus[D];
	
	for (int d = 0; d < D; ++d) {
		d_minus[d] = (log(x[d] / Strike) + (r - delta - sigma * sigma / 2.) * (T - t)) / sigma 	/ sqrt(T - t + 0.000001);
		d_plus[d] = d_minus[d] + sigma * sqrt(T - t); 
	}

	double erg = 0;
	for (int l = 0; l < D; ++l) {
		double integralSumme = 0;
		double dz = dt;
		
		for (double z = -5; z < d_plus[l]; z += dz) {
			double df = exp(-0.5 * z * z);
			for (int l_Strich = 0; l_Strich < D; ++l_Strich)
				if (l_Strich != l)
					df *= cnd(log(x[l] / x[l_Strich]) / sigma / sqrt(T - t + 0.000001) - z + sigma * sqrt(T - t));
				
			integralSumme += df * dz;
		}

		erg += x[l] * exp(-delta * (T - t)) / sqrt(2 * 3.141592654)	* integralSumme;
	}
	
	double prod = 1;
	for (int l = 0; l < D; ++l)
		prod *= (1 - cnd(d_minus[l]));
	erg += -Strike * exp(-r * (T - t)) + Strike * exp(-r * (T - t)) * prod;
	
	return exp(-r * t) * (erg);
}


EuroBewerter::EuroBewerter() {
}

EuroBewerter::~EuroBewerter() {	
}

}
