#include "AmericanOption.h"

#include <stdlib.h>
#include <string.h>

using namespace std;

AmericanOption::AmericanOption() {
	Daten();
	dt = T / (double) (N - 1);
}

AmericanOption::~AmericanOption() {
	delete[] X0; // Spot
	delete[] sigma; //Volatility
}

double AmericanOption::payoff(double* x, int time) {
	if (option == MIN_PUT)
		return std::max(Strike - Min(x, D), 0.) * exp(-r * dt * (double) (time)); //Min Put
	if (option == MAX_CALL)
		return std::max(Max(x, D) - Strike, 0.) * exp(-r * dt * (double) (time)); //Max Call

	printf("ERROR, option unknown!\n");
	exit(0);
	return -1;
}

void AmericanOption::Pfadgenerieren(double** X, int start, double* S,
		RNG* generator) {
	double** wdiff = DoubleFeld(N, D);
	for (int n = 0; n < N; ++n)
		for (int d = 0; d < D; ++d)
			wdiff[n][d] = sqrt(dt) * generator->nextGaussian();
	Pfadgenerieren(X, wdiff, start, S);
	deleteDoubleFeld(wdiff, N, D);
}

void AmericanOption::Pfadgenerieren(double** X, double** wdiff, int start,
		double * S) {
	for (int d = 0; d < D; ++d)
		X[start][d] = S[d];
	for (int d = 0; d < D; ++d) {
		for (int n = start + 1; n < N; ++n) {
			X[n][d] = X[n - 1][d]
					* exp(
							(((r - delta) - 0.5 * sigma[d] * sigma[d]) * dt
									+ sigma[d] * wdiff[n][d]));
		}
	}
}
