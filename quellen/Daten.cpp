#include "AmericanOption.h"

using namespace std;

void AmericanOption::Daten() {
	D = 2;

	X0 = DoubleFeld(D); 
	sigma = DoubleFeld(D);
		 
	for (int j = 0; j < D; ++j) {
		X0[j] = 90.;
		sigma[j] = 0.2;
	}
	
	delta = 0.1;
	Strike = 100.;
	r = 0.05;
	T = 3;
	N = 10;
}

