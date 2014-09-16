#include <stdio.h>
#include "AmericanOption.h"
#include <cstring>
#include <string.h>
#include <iostream>
#include <algorithm>
#include "Linear_Regression.h"
#include <time.h>
#include <assert.h>
#include <iostream>   // std::cout
#include <string>     // std::string, std::to_string
#include <time.h>
#include "Hilfsmittel.h"
#include "AmericanOption.h"
using namespace std;

void AmericanOption::trainingpaths_erstellen(int l) {
	RNG generator;
	//generator.setSeed(7);
	double** wdiff = DoubleFeld(N, D);
	for (int m = 0; m < Mtraining(l); ++m) {
		if (l == 0 || m >= Mtraining(l-1)) {
			for (int d = 0; d < D; ++d)
				for (int n = 0; n < N; ++n)
					if (m % 2 == 0)
						wdiff[n][d] = generator.nextGaussian() * sqrt(dt);
					else
						wdiff[n][d] *= -1.;
			Pfadgenerieren(X[l][m], wdiff, 0, X0);
		} else
			for (int d = 0; d < D; ++d)
				for (int n = 0; n < N; ++n)
					X[l][m][n][d] = X[l - 1][m][n][d];
	}
	deleteDoubleFeld(wdiff, N, D);
}

void AmericanOption::ErgebnisseAusgeben(int l) {
	char buffer[50];

	sprintf(buffer, "mittelwert%d.txt", l);
	ofstream File(buffer, ios::out | ios::app);
	if (File.is_open())
		File << mittelwert(Levelergs[l]) << endl;

	sprintf(buffer, "varianz%d.txt", l);
	ofstream File2(buffer, ios::out | ios::app);
	if (File2.is_open())
		File2 << varianz(Levelergs[l]) << endl;

}

void AmericanOption::printInfo() {
	double sum = 0;
	for (int l = 0; l < L; ++l) {
		printf("Level %d:\t %.5lf (%.5lf) \t %d Pfade", l,
				mittelwert(Levelergs[l]), varianz(Levelergs[l]),
//				RegressionV(Levelergs[l],cv[l],EB.european_MaxCall_ND(X0,D,0,T,Strike,r,delta,sigma[0],0.001)),
				(int) Levelergs[l].size());
		if (n != NULL)
			printf(" n=%d ", n[l]);
		printf("\n");
		sum += mittelwert(Levelergs[l]);
	}
	printf("sum=%f\n", sum);
	printf("\n\n");
}

void AmericanOption::addLevelPath(int l) {
	static RNG generator;
	double** XX = NULL;
	double ** wdiff = NULL;

	wdiff = DoubleFeld(N, D);
	XX = DoubleFeld(N, D);

	for (int d = 0; d < D; ++d)
		for (int n = 0; n < N; ++n)
			wdiff[n][d] = generator.nextGaussian() * sqrt(dt);
	Pfadgenerieren(XX, wdiff, 0, X0);
	double e1 = Pfad(XX, l);

//	for (int d = 0; d < D; ++d)
//		for (int n = 0; n < N; ++n)
//			wdiff[n][d] *= -1.;
//	Pfadgenerieren(XX, wdiff, 0, X0);
//	double e2 = Pfad(XX, l);

	Levelergs[l].push_back(e1);

//	Levelergs[l].push_back(0.5 * (e1 + e2));

	deleteDoubleFeld(XX, N, D);
	deleteDoubleFeld(wdiff, N, D);
}

double AmericanOption::Pfad(double ** XX, int l) {
	double erg = 0;
//int stopp=N-1;

	for (int n = 0; n < N; ++n) {
		if (payoff(XX[n], n) > 0 || n == N - 1)
			if (payoff(XX[n], n) >= C_estimate_Mesh(XX[n], n, l)) {
				erg += payoff(XX[n], n);
				erg -=
						1.2
								* (EB.european_MaxCall_ND(XX[n], D,
										(double) (n) * dt, T, Strike, r, delta,
										sigma[0], 0.005) - 6.65509);
//				erg -= 1.2
//						* (EB.european_MaxCall_ND(XX[n], D, (double) (n) * dt,
//								T, Strike, r, delta, sigma[0], 0.01)
//								- EB.european_MaxCall_ND(X0, D, 0, T, Strike, r,
//										delta, sigma[0], 0.01));

				break;
			}
	}

	if (l > 0)
		for (int n = 0; n < N; ++n) {
			if (payoff(XX[n], n) > 0 || n == N - 1)
				if (payoff(XX[n], n) >= C_estimate_Mesh(XX[n], n, l - 1)) {
					erg -= payoff(XX[n], n);
					erg += 1.2
							* (EB.european_MaxCall_ND(XX[n], D,
									(double) (n) * dt, T, Strike, r, delta,
									sigma[0], 0.005) - 6.65509);
//					erg -= 1.2
//							* (EB.european_MaxCall_ND(XX[n], D,
//									(double) (n) * dt, T, Strike, r, delta,
//									sigma[0], 0.01)
//									- EB.european_MaxCall_ND(X0, D, 0, T,
//											Strike, r, delta, sigma[0], 0.01));

					break;
				}
		}

	return erg;

}

void AmericanOption::simplified() {
	printf("simplified multilevel\n");

	Daten();

	double eps = 0.05;
	int maxL = 15;
	n = new int[maxL];

	dt = T / (double) (N - 1);
	X = DoubleFeld(maxL, Mtraining(maxL - 1), N, D);
	V = DoubleFeld(maxL, Mtraining(maxL - 1), N);
	weights = DoubleFeld(maxL, N, Mtraining(maxL - 1));
	weight_sum = DoubleFeld(maxL, N, Mtraining(maxL - 1));

	Levelergs = new vector<double> [maxL];

	L = 1;

	bool done = false;

	double kappa2 = 1.;

	while (!done) {

		trainingpaths_erstellen(L - 1);
		weights_erstellen(L - 1);
		trainingpaths_regression(L - 1);

		for (int ii = 0; ii < 100; ++ii) {
			addLevelPath(L - 1);
			if (ii % 1000 == 0)
				printInfo();
		}

		double vorder = 0;
		for (int l = 0; l < L; ++l)
			vorder += sqrt(
					pow((double) Mtraining(l), kappa2) * varianz(Levelergs[l]));

		for (int ll = 0; ll < L; ++ll) {
			double hinter = sqrt(
					pow((double) Mtraining(ll), -kappa2)
							* varianz(Levelergs[ll]));
			n[ll] = ceil(3. / eps / eps * vorder * hinter);
			//	printf("vorder= %.30lf, hinter=%.30lf\n", vorder, hinter);
			printf("\n\nn[%d]=%d\n", L - 1, n[L - 1]);
		}
//		sleep(2);
		done = true;

		if (L >= 3) {
			double Pl = 0;
			for (int ll = 0; ll < L; ++ll) {
				for (int i = Levelergs[ll].size(); i < n[ll]; ++i) {
					addLevelPath(ll);
					if (i % 1000 == 0 || i == n[ll] - 1)
						printInfo();
				}
				Pl += mittelwert(Levelergs[ll]);
			}

//			if (mittelwert(Levelergs[L - 1]) / Pl > eps / sqrt(3.)) { //old criterion
			if (max((mittelwert(Levelergs[L - 1])),
					0.5 * (mittelwert(Levelergs[L - 2]))) > eps / sqrt(3.)) { //new criterion
				L += 1;
				done = false;
			}

		}
		if (L < 3) {
			L += 1;
			done = false;
		}
	}

	double summe = 0;
	for (int l = 0; l < L; ++l)
		summe += mittelwert(Levelergs[l]);

	ofstream File("simpergs.txt", ios::out | ios::app);
	if (File.is_open())
		File << summe << endl;

	ofstream FileLN("numberOfLevels.txt", ios::out | ios::app);
	if (FileLN.is_open())
		FileLN << L << endl;

	long double comp_training1 = 0;
	long double comp_training2 = 0;
	long double comp_testing = 0;

	for (int l = 0; l < L; ++l) {
		comp_training1 += pow(Mtraining(l), kappa2 + 1);
		if (l > 0)
			comp_training1 += pow(Mtraining(l - 1), kappa2 + 1);
		comp_training2 += pow(Mtraining(l), kappa2 + 1);
		comp_testing += n[l] * pow(Mtraining(l), kappa2);
		if (l > 0)
			comp_testing += n[l] * pow(Mtraining(l - 1), kappa2);
	}

//	ofstream FileC1("comp_training.txt", ios::out | ios::app);
//	if (FileC1.is_open())
//		FileC1 << comp_training1 << endl;

	ofstream FileC2("comp_training_trick.txt", ios::out | ios::app);
	if (FileC2.is_open())
		FileC2 << comp_training2/1000000. << endl;

	ofstream FileC3("comp_testing.txt", ios::out | ios::app);
	if (FileC3.is_open())
		FileC3 << comp_testing/1000000. << endl;

}

void AmericanOption::run() {
	Daten();

	L = 6;
	dt = T / (double) (N - 1);
	X = DoubleFeld(L, Mtraining(L - 1), N, D);
	V = DoubleFeld(L, Mtraining(L - 1), N);
	weights = DoubleFeld(L, N, Mtraining(L - 1));
	weight_sum = DoubleFeld(L, N, Mtraining(L - 1));

	n = NULL;
	Levelergs = new vector<double> [L];

	for (int l = 0; l < L; ++l) {
		printf("Level %d: %d training paths \n", l, Mtraining(l));
		trainingpaths_erstellen(l);
		weights_erstellen(l);
		trainingpaths_regression(l); // Mesh
	}

	for (int kk = 0; kk < 10; ++kk) {
		for (int i = 0; i < 100; ++i)
			for (int l = 0; l < L; ++l)
				addLevelPath(l);
		printInfo();
	}

	for (int l = 0; l < L; ++l)
		ErgebnisseAusgeben(l);
}

int main(int argc, char* args[]) {
	AmericanOption AMO;
//long double t=435343453454345345;
//cout.precision(100);
//
//ofstream FileC2("trick.txt", ios::out | ios::app);
//	if (FileC2.is_open())
//		FileC2 <<  << endl;
//
//
//	exit(0);

	if (argc > 1) {
		string arg = args[1];
		if (!arg.compare("-simp"))
			AMO.simplified();
	} else
		AMO.run();
}

