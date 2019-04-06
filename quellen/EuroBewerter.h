#ifndef EUROBEWERTER_H_
#define EUROBEWERTER_H_


#include <functional>

namespace std {

class EuroBewerter {
public:
	EuroBewerter();
	virtual ~EuroBewerter();

	double call(double t, double T, double X0, double Strike, double r, double delta, double sigma);
	double put(double t, double T, double X0, double Strike, double r, double delta, double sigma);
	double exchange_option(double x, double y, double t, double T, double r, double delta, double sigma);
	double european_MaxCall_ND(double* x, int D, double t, double T,double Strike, double r, double delta, double sigma, double dt);
};
} 
#endif 
