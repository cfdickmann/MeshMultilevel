#include "AmericanOption.h"
//#include "../alglib/solvers.h"
#include "Hilfsmittel.h"
#include "Linear_Regression.h"

using namespace std;

int AmericanOption::Mtraining(int l) {
	return 20 * pow(2, l);
}

void AmericanOption::trainingpaths_regression(int l) {
	for (int n = N - 1; n > 0; --n) {
		for (int m = 0; m < Mtraining(l); ++m) {
			if (m % 10 == 0)
				printf("Mesh Training Schritt %d %.0lf \%% \r", n,
						(double) m / (double) Mtraining(l) * 100.);
			cout.flush();
			double est = C_estimate_Mesh(X[l][m][n], n, l);
			V[l][m][n] = max(est, payoff(X[l][m][n], n));
		}
	}
}

void AmericanOption::weights_erstellen(int ll) {
	for (int n = 0; n < N - 1; ++n)
		for (int l = 0; l < Mtraining(ll); ++l)
			weights[ll][n][l] = kernelD(X0, X[ll][l][n + 1],
					dt * (double) (n + 1));

	for (int n = 0; n < N - 1; ++n)
		for (int l = 0; l < Mtraining(ll); ++l) {
			weight_sum[ll][n][l] = 0;
			for (int m = 0; m < Mtraining(ll); ++m)
				weight_sum[ll][n][l] += kernelD(X[ll][m][n], X[ll][l][n + 1],
						dt);
		}
}

double phi(double x) {
	static double FF = 1. / sqrt(2. * 3.141592653);
	return FF * exp(-(x * x) / 2.);
}

double AmericanOption::kernelD(double* von, double* nach, double dt) {
	double produkt = 1.;
	for (int d = 0; d < D; ++d) {
		double ratio = nach[d] / von[d];
		double klammer = (log(ratio)
				- (r - delta - 0.5 * sigma[d] * sigma[d]) * dt)
				/ (sigma[d] * sqrt(dt));
		produkt *= phi(klammer) / (sigma[d] * sqrt(dt) * nach[d]);
	}
	return produkt;
}

double AmericanOption::C_estimate_Mesh(double* x, int Etime, int l) {
	if (Etime == N - 1)
		return 0;

	double int_dt = 0.005;

	double v = EB.european_MaxCall_ND(x, D, (double) (Etime) * (dt), //gutes
	(double) (Etime + 1) * dt, Strike, r, delta, sigma[0], int_dt);

	vector<double> weightsss;
	vector<double> y;
	vector<double> cv;
	for (int k = 0; k < Mtraining(l); ++k) {
		double K = kernelD(x, X[l][k][Etime + 1], dt) / weight_sum[l][Etime][k];
		double CV = exp(-r * (double) (Etime + 1) * dt)
				* std::max(Max(X[l][k][Etime + 1], D) - Strike, 0.);	// gutes
		weightsss.push_back(K / (double) (Mtraining(l) * Mtraining(l)));
		y.push_back(V[l][k][Etime + 1]);
		cv.push_back(CV);
	}

	if (int_dt == 0)
		return mittelwert(y); //gutes
	return RegressionV(y, cv, weightsss, v); //gutes
}
