#include "EuroBewerter.h"
#include "math.h"
#include <stdio.h>
#include <functional>
#include <assert.h>
#include <stdlib.h>
#define pi 3.141592654

using namespace std;

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

    double xGL[25]= {
        1,
0.997005175362643586236,
0.989972220058020802022,
0.9789519106467198092673,
0.9639894801099686121339,
0.9451453286808814137682,
0.922495378502752746657,
0.8961308490768326327616,
0.8661579081773533109677,
0.832697249796457952787,
0.795883610117728469387,
0.7558652257940530390604,
0.7128032374236458598544,
0.6668710408790448158956,
0.6182535892015291708338,
0.5671466479156704257289,
0.5137560067830487440676,
0.4582966511792997592748,
0.4009918964363984808992,
0.3420724886388038436291,
0.281775675495909512023,
0.2203442510330564760657,
0.1580255779484029029865,
0.0950705915726566477616,
0.0317327894426237129119,
    };
    
    double wGL [25]={
        8.163265306122448979592E-4,
0.0050273455497416577264,
0.0090329866868996020993,
0.013000139812211760657,
0.0169146554838637660232,
0.0207609894970159934393,
0.024523701718169494315,
0.0281876554510602491663,
0.0317381026065706387509,
0.03516074908934724962954,
0.0384418142109977345625,
0.041568086852368492733,
0.04452697893447274223,
0.04730657622923672285351,
0.04989568639334337604645,
0.052283884066433463416,
0.054461552866930615293,
0.0564199241232101410276,
0.058151112187474100161,
0.0596481461918334541467,
0.060904998119632791003,
0.0619166070794861790183,
0.06267889968456801074686,
0.0631888064552666332142,
0.0634442741792528868443
    };
    


static double IntegrationGaussLobatto(std::function<double(double)> func)
{    
    double integral =0.0;
    
    for(int i=0; i<25; ++i)
        integral+= func(-xGL[i]) * wGL[i];
    
    for(int i=0; i<25; ++i)
        integral+= func(xGL[i]) * wGL[i];
    
    return integral;
}

static double IntegrationGaussLobattoInterval(std::function<double(double)> func, double a, double b)
{
     std::function<double(double)> f =[&](double x)
    {
        return func(a + (b-a) * (x + 1.0) / 2.0);
    };
    
    return IntegrationGaussLobatto(f)*(b-a)/2.0;        
}

static double testInt()
{
    std::function<double(double)> f =[&](double x)
    {
        return exp(-x*x);
    };
    
    double erg=IntegrationGaussLobatto(f);

    assert( abs(erg-1.49365)<0.00001);
}

static double testInt2()
{
    std::function<double(double)> f =[&](double x)
    {
        return exp(-x*x);
    };
    
    double erg=IntegrationGaussLobattoInterval(f,-2,3);

    assert( abs(erg-1.76829)<0.00001);
}

double EuroBewerter::european_MaxCall_ND(double* x, int D, double t, double T,
		double Strike, double r, double delta, double sigma, double dt) {
	
	if(D==1)
		return call(t,T,x[0],Strike,r,delta,sigma);
	
	if(dt==0)
		return 0;    
    
    //testInt();
    //testInt2();
    
	double d_minus[D];
	double d_plus[D];
	
	for (int d = 0; d < D; ++d) {
		d_minus[d] = (log(x[d] / Strike) + (r - delta - sigma * sigma / 2.) * (T - t)) / sigma 	/ sqrt(T - t + 0.000001);
		d_plus[d] = d_minus[d] + sigma * sqrt(T - t); 
	}

	double erg = 0;
	for (int l = 0; l < D; ++l) {	
		
		 std::function<double(double)> fu =[&](double z)
        {
            double df = exp(-0.5 * z * z);
                for (int l_Strich = 0; l_Strich < D; ++l_Strich)
                    if (l_Strich != l)
                        df *= cnd(log(x[l] / x[l_Strich]) / sigma / sqrt(T - t + 0.000001) - z + sigma * sqrt(T - t));
                
                    return df;
        };
		
		double integralSumme=IntegrationGaussLobattoInterval(fu, -5.0, d_plus[l]);

		erg += x[l] * exp(-delta * (T - t)) / sqrt(2 * 3.141592654)	* integralSumme;
	}
	
	double prod = 1;
	for (int l = 0; l < D; ++l)
		prod *= (1 - cnd(d_minus[l]));
	erg += -Strike * exp(-r * (T - t)) + Strike * exp(-r * (T - t)) * prod;
	
	return exp(-r * t) * (erg);
}

/*
double EuroBewerter::european_MaxCall_ND(double* x, int D, double t, double T,
		double Strike, double r, double delta, double sigma, double dt) {
	
	if(D==1)
		return call(t,T,x[0],Strike,r,delta,sigma);
	
	if(dt==0)
		return 0;    
    
    testInt();
    
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
}*/


EuroBewerter::EuroBewerter() {
}

EuroBewerter::~EuroBewerter() {	
}

}
