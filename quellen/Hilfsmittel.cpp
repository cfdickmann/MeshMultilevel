#include "Hilfsmittel.h"
#include <stdio.h>
#include <iostream>
#include <string.h>
#include <math.h>
#include <iostream>
#include <fstream>
#include <stdlib.h>
#include <cstring>
#include <math.h>
#include <iostream>
#include <algorithm>
#include <functional>
#include <vector>
#include <stdlib.h>
//#include "../alglib/ap.h"
//#include "../alglib/linalg.h"
//#include "../alglib/solvers.h"
#include "AmericanOption.h"
//using namespace alglib;
using namespace std;
//
//double varianz(vector<double> vec){
//	if(vec.size()==0)return 0;
//	double erg=0;
//	double ergQ;
//	for(unsigned i=0;i<vec.size();++i)
//	{
//		erg+=vec.at(i);
//		ergQ+=pow(vec.at(i),2);
//	}
//	return (ergQ/(double)vec.size())-pow(erg/(double)vec.size(),2);
//}

double varianz(vector<double> vec) {
	if (vec.size() == 0)
		return 0;
	double erg = 0;
	double m = mittelwert(vec);
	for (unsigned i = 0; i < vec.size(); ++i) {
		erg += pow(vec.at(i) - m, 2);
	}
	return (erg) / ((double) vec.size() - 1.);
}

double mittelwert(vector<double> vec) {
	if (vec.size() == 0)
		return 0;
	double erg = 0;
	for (unsigned i = 0; i < vec.size(); ++i)
		erg += vec.at(i);
	return erg / vec.size();
}

void MatrixAusgeben(double** a, int D) {
	printf("\n");
	for (int d = 0; d < D; ++d) {
		for (int f = 0; f < D; ++f)
			printf("%.10lf, ", a[d][f]);
		printf("\n");
	}
}

double** MatrixMultiplizieren(double** a, double** b, int D) {
	double** erg = DoubleFeld(D, D);
	for (int d = 0; d < D; ++d)
		for (int f = 0; f < D; ++f) {
			double summe = 0;
			for (int k = 0; k < D; ++k)
				summe += a[d][k] * b[k][f];
			erg[d][f] = summe;
		}
	return erg;
}

double betrag(double x) {
	return x < 0 ? -x : x;
}

int binary(int number, int digit) {
	int rest = number % (int) pow(2, digit + 1);
	return rest / (int) pow(2, digit);
}

double* alphasLaden(int K) {
	double* alpha = (double*) malloc(sizeof(double) * (K + 1));
	string line;
	ifstream myfile("alphas.dat");
	int k = 0;
	if (myfile.is_open()) {
		while ((!myfile.eof()) && k < K) {
			//printf("gelesen\n");
			getline(myfile, line);
			char * buffer = new char[line.length()];
			strcpy(buffer, line.c_str());
			alpha[k] = (double) (atof(buffer));
			k++;
		}
		myfile.close();
		int temp = k - 1;
		for (int i = k; i < K; ++i)
			alpha[i] = alpha[(i - 1) % temp + 1];
		return alpha;
	} else
		printf("Kann Datei nicht oeffnen\n");
	return alpha;
}

void werteSchreiben(double* w, int K, int N) {
	fstream f;
	f.open("werte.dat", ios::out);
	f << w[0];
	for (int i = 1; i < K; ++i) {
		if (i % N == 0)
			f << endl;
		f << " " << w[i];
	}
	f.close();
}

void alphasSchreiben(double* alpha, int K) {
	fstream f;
	f.open("alphas.dat", ios::out);
	f << floor(10000. * alpha[0]) / 10000.;
	for (int i = 1; i < K; ++i)
		f << endl << floor(10000. * alpha[i]) / 10000.;
	f.close();
}

int argMin(double* v, int l) {
	double minimum = v[0];
	int minimum_j = 0;
	for (int j = 1; j < l; ++j)
		if (v[j] < minimum) {
			minimum_j = j;
			minimum = v[j];
		}
	return minimum_j;
}

void ErgebnisAnhaengen(double d) {
	ofstream File("ergebnisse2.txt", ios::out | ios::app);
	if (File.is_open())
		File << d << endl;
}

void ErgebnisAnhaengenML(double d) {
	ofstream File("ergebnisseML2.txt", ios::out | ios::app);
	if (File.is_open())
		File << d << endl;
}

void ErgebnisAnhaengen(vector<double> d, char* filename) {
	ofstream File(filename, ios::out | ios::app);
	if (File.is_open()) {
		for (unsigned i = 0; i < d.size(); ++i)
			File << d.at(i) << endl;
	}
}

void ErgebnisAnhaengen(double d, char* filename) {
	ofstream File(filename, ios::out | ios::app);
	if (File.is_open())
		File << d << endl;
}

void BubbleSort(double* werte, int* index, int l) {
	bool geaendert = true;
	while (geaendert) {
		geaendert = false;
		for (int i = 1; i < l; ++i) {
			if (werte[index[i]] > werte[index[i - 1]]) {
				int temp = index[i - 1];
				index[i - 1] = index[i];
				index[i] = temp;
				geaendert = true;
			}
		}
	}
}

int argMax(double* v, int l) {
	double maximum = v[0];
	int maximum_j = 0;
	for (int j = 1; j < l; ++j)
		if (v[j] > maximum) {
			maximum_j = j;
			maximum = v[j];
		}
	return maximum_j;
}

int argZweiter(double *v, int l) {
	double maximum = -10000000;
	double max = argMax(v, l);

	int argZwe = 0;
	for (int j = 0; j < l; ++j)
		if (v[j] > maximum && j != max) {
			argZwe = j;
			maximum = v[j];
		}
	return argZwe;
}

int argDritter(double *v, int l) {
	double maximum = -10000000;
	double max = argMax(v, l);
	double zwe = argZweiter(v, l);
	int argZwe = 0;
	for (int j = 0; j < l; ++j)
		if (v[j] > maximum && j != max && j != zwe) {
			argZwe = j;
			maximum = v[j];
		}
	return argZwe;
}

void InPipeSchreiben(int* pipe, double wert) {
//	close(pipe[0]);
//	char string[20];
//	sprintf(string, "%f", wert);
//	write(pipe[1], string, (strlen(string) + 1));
}
//
//double AusPipeLesen(int* pipe) {
////	close(pipe[1]);
////	char readbuffer[200];
////	read(pipe[0], readbuffer, sizeof(readbuffer));
////	return atof(readbuffer);
//}

double Max(double* v, int l) {
	return v[argMax(v, l)];
}

double Min(double* v, int l) {
	return v[argMin(v, l)];
}

////double qnorm(double p) {
//	/** * @(#)qnorm.js * * Copyright (c) 2000 by Sundar Dorai-Raj
//	 * * @author Sundar Dorai-Raj
//	 * * Email: sdoraira@vt.edu
//	 * * This program is free software; you can redistribute it and/or
//	 * * modify it under the terms of the GNU General Public License
//	 * * as published by the Free Software Foundation; either version 2
//	 * * of the License, or (at your option) any later version,
//	 * * provided that any use properly credits the author.
//	 * * This program is distributed in the hope that it will be useful,
//	 * * but WITHOUT ANY WARRANTY; without even the implied warranty of
//	 * * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the
//	 * * GNU General Public License for more details at http://www.gnu.org * * */
//	// ALGORITHM AS 111, APPL.STATIST., VOL.26, 118-121, 1977.
//	// Computes z=invNorm(p)
//	//	double split = 0.42;
//	//	double a0 = 2.50662823884, a1 = -18.61500062529, a2 = 41.39119773534;
//	//	double a3 = -25.44106049637, b1 = -8.47351093090, b2 = 23.08336743743;
//	//	double b3 = -21.06224101826, b4 = 3.13082909833;
//	//	double c0 = -2.78718931138, c1 = -2.29796479134, c2 = 4.85014127135;
//	//	double c3 = 2.32121276858, d1 = 3.54388924762, d2 = 1.63706781897;
//	//	double q = p - 0.5;double rr = 0;double ppnd = 0;
//	//	if (fabs(q) <= split) {
//	//		rr = q * q;
//	//		ppnd = q * (((a3 * rr + a2) * rr + a1) * rr + a0) / ((((b4 * rr + b3) * rr + b2) * rr + b1) * rr + 1);
//	//	} else {
//	//		rr = p;
//	//		if (q > 0)
//	//			rr = 1 - p;
//	//		if (rr > 0) {
//	//			rr = sqrt(-log(rr));
//	//			ppnd = (((c3 * rr + c2) * rr + c1) * rr + c0) / ((d2 * rr + d1) * rr + 1);
//	//			if (q < 0)
//	//				ppnd = -ppnd;
//	//		} else
//	//			ppnd = 0;
//	//	}
//	//	return (ppnd);
//
//	/** * @(#)qnorm.js * * Copyright (c) 2000 by Sundar Dorai-Raj
//	 * * @author Sundar Dorai-Raj
//	 * * Email: sdoraira@vt.edu
//	 * * This program is free software; you can redistribute it and/or
//	 * * modify it under the terms of the GNU General Public License
//	 * * as published by the Free Software Foundation; either version 2
//	 * * of the License, or (at your option) any later version,
//	 * * provided that any use properly credits the author.
//	 * * This program is distributed in the hope that it will be useful,
//	 * * but WITHOUT ANY WARRANTY; without even the implied warranty of
//	 * * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the
//	 * * GNU General Public License for more details at http://www.gnu.org * * */
//	// ALGORITHM AS 111, APPL.STATIST., VOL.26, 118-121, 1977.
//	// Computes z=invNorm(p)
//	double split = 0.42;
//	double a0 = 2.50662823884, a1 = -18.61500062529, a2 = 41.39119773534;
//	double a3 = -25.44106049637, b1 = -8.47351093090, b2 = 23.08336743743;
//	double b3 = -21.06224101826, b4 = 3.13082909833;
//	double c0 = -2.78718931138, c1 = -2.29796479134, c2 = 4.85014127135;
//	double c3 = 2.32121276858, d1 = 3.54388924762, d2 = 1.63706781897;
//	double q = p - 0.5;double rr = 0;double ppnd = 0;
//	if (abs(q) <= split) {
//		rr = q * q;
//		ppnd = q * (((a3 * rr + a2) * rr + a1) * rr + a0) / ((((b4 * rr + b3) * rr + b2) * rr + b1) * rr + 1);
//	} else {
//		rr = p;
//		if (q > 0)
//			rr = 1 - p;
//		if (rr > 0) {
//			rr = sqrt(-log(rr));
//			ppnd = (((c3 * rr + c2) * rr + c1) * rr + c0) / ((d2 * rr + d1) * rr + 1);
//			if (q < 0)
//				ppnd = -ppnd;
//		} else
//			ppnd = 0;
//	}
//	return (ppnd);
//}

double max(double x, double y) {
	return x < y ? y : x;
}

void ausgeben(double* x, int j) {
	printf("[");
	for (int k = 0; k < j - 1; ++k)
		printf("%.4lf, ", x[k]);
	printf("%.4lf", x[j - 1]);
	printf("]\n");
}

double CumulativeNormalDistribution(double x) {
	int neg = (x < 0.) ? 1 : 0;
	if (neg == 1)
		x *= -1.;
	double k = (1. / (1. + 0.2316419 * x));
	double y = ((((1.330274429 * k - 1.821255978) * k + 1.781477937) * k
			- 0.356563782) * k + 0.319381530) * k;
	y = 1.0 - 0.398942280401 * exp(-0.5 * x * x) * y;
	return (1. - neg) * y + neg * (1. - y);
}

int * IntFeld(int m) {
	int* erg = new int[m];
	return erg;
}

int ** IntFeld(int m, int n) {
	int** erg = new int*[m];
	for (int i = 0; i < m; ++i)
		erg[i] = IntFeld(n);
	return erg;
}

int *** IntFeld(int m, int n, int o) {
	int*** erg = new int**[m];
	for (int i = 0; i < m; ++i)
		erg[i] = IntFeld(n, o);
	return erg;
}

double * DoubleFeld(int m) {
	double* erg = new double[m];
	return erg;
}

double ** DoubleFeld(int m, int n) {
	double** erg = new double*[m];
	for (int i = 0; i < m; ++i)
		erg[i] = new double[n];
	return erg;
}

double *** DoubleFeld(int m, int n, int o) {
	double*** erg = new double**[m];
	for (int i = 0; i < m; ++i)
		erg[i] = DoubleFeld(n, o);
	return erg;
}

double **** DoubleFeld(int m, int n, int o, int p) {
	double **** erg = new double ***[m];
	for (int i = 0; i < m; ++i)
		erg[i] = DoubleFeld(n, o, p);
	return erg;
}

double ***** DoubleFeld(int m, int n, int o, int p, int q) {
	double ***** erg = new double****[m];
	for (int i = 0; i < m; ++i)
		erg[i] = DoubleFeld(n, o, p, q);
	return erg;
}

void deleteDoubleFeld(double * D, int m) {
	delete[] D;
}

void deleteDoubleFeld(double ** D, int m, int n) {
	for (int i = 0; i < m; ++i)
		deleteDoubleFeld(D[i], n);
	delete[] D;
}

void deleteDoubleFeld(double *** D, int m, int n, int o) {
	for (int i = 0; i < m; ++i)
		deleteDoubleFeld(D[i], n, o);
	delete[] D;
}

void deleteDoubleFeld(double **** D, int m, int n, int o, int p) {
	for (int i = 0; i < m; ++i)
		deleteDoubleFeld(D[i], n, o, p);
	delete[] D;
}

void deleteDoubleFeld(double ***** D, int m, int n, int o, int p, int q) {
	for (int i = 0; i < m; ++i)
		deleteDoubleFeld(D[i], n, o, p, q);
	delete[] D;
}

void deleteIntFeld(int * D, int m) {
	delete[] D;
}

void deleteIntFeld(int ** D, int m, int n) {
	for (int i = 0; i < m; ++i)
		deleteIntFeld(D[i], n);
	delete[] D;
}

void deleteIntFeld(int *** D, int m, int n, int o) {
	for (int i = 0; i < m; ++i)
		deleteIntFeld(D[i], n, o);
	delete[] D;
}

int* array_machen(int z) {
	int* erg = (int*) malloc(sizeof(int));
	erg[0] = z;
	return erg;
}

//int* AmericanOption::pivot(double** A) {
//	//	int nn = Mphi;
//	int* pivot = (int*)malloc(sizeof(int)*Mphi);
//	for (int j = 0; j < Mphi - 1; j++) {
//		double max = fabs(A[j][j]);
//		int imax = j;
//		for (int i = j + 1; i < Mphi; i++)
//			if (fabs(A[i][j]) > max) {
//				max = fabs(A[i][j]);
//				imax = i;
//			}
//		double* h =DoubleFeld(Mphi);
//		for(int i=0;i<Mphi;++i)h[i]=A[j][i];
//		A[j] = A[imax];
//		A[imax] = h;
//		pivot[j] = imax;
//		for (int i = j + 1; i < Mphi; i++) {
//			double f = -A[i][j] / A[j][j];
//			for (int k = j + 1; k < Mphi; k++)
//				A[i][k] += f * A[j][k];
//			A[i][j] = -f;
//		}
//	}
//	return pivot;
//}

//int* AmericanOption::pivot(double** A) {
//	//	int nn = Mphi;
//	int* pivot = new int[Mphi];
//	for (int j = 0; j < Mphi - 1; j++) {
//		double max = fabs(A[j][j]);
//		int imax = j;
//		for (int i = j + 1; i < Mphi; i++)
//			if (fabs(A[i][j]) > max) {
//				max = fabs(A[i][j]);
//				imax = i;
//			}
//		double* h = DoubleFeld(Mphi);
//		for (int i = 0; i < Mphi; ++i)
//			h[i] = A[j][i];
//		A[j] = A[imax];
//		A[imax] = h;
//		pivot[j] = imax;
//		for (int i = j + 1; i < Mphi; i++) {
//			double f = -A[i][j] / A[j][j];
//			for (int k = j + 1; k < Mphi; k++)
//				A[i][k] += f * A[j][k];
//			A[i][j] = -f;
//		}
//	}
//	return pivot;
//}
//
//double* AmericanOption::LGSloesen(double** AA, double* bb, int Mphi) {
//	//return LGS_mit_alglib_loesen(AA,bb,Mphi);
//
//	double** A = new double*[Mphi];
//	for (int m = 0; m < Mphi; ++m)
//		A[m] = new double[Mphi];
//
//	for (int m = 0; m < Mphi; ++m)
//		for (int n = 0; n < Mphi; ++n)
//			A[m][n] = AA[m][n];
//
//	double* b = new double[Mphi];
//	for (int m = 0; m < Mphi; ++m)
//		b[m] = bb[m];
//	double** B = DoubleFeld(Mphi, Mphi);
//	for (int i = 0; i < Mphi; ++i)
//		for (int j = 0; j < Mphi; ++j)
//			B[i][j] = A[i][j];
//
//	//	double** A=AA;
//	//	double* b=bb;
//	//	double** B=A;
//	//	// loest das LGS Ax = b nach x auf
//
//	double* x = DoubleFeld(Mphi);
//	for (int i = 0; i < Mphi; ++i)
//		x[i] = b[i];
//
//	int* piv = pivot(B);
//	int nn = Mphi;
//	for (int i = 0; i < nn - 1; i++) {
//		double h = b[piv[i]];
//		b[piv[i]] = b[i];
//		b[i] = h;
//	}
//	for (int j = 0; j < nn; j++) {
//		x[j] = b[j];
//		for (int i = 0; i < j; i++)
//			x[j] -= B[j][i] * x[i];
//	}
//	for (int j = nn - 1; j >= 0; j--) {
//		for (int k = j + 1; k < nn; k++)
//			x[j] -= B[j][k] * x[k];
//		x[j] /= B[j][j];
//	}
//	deleteIntFeld(piv, Mphi);
//	deleteDoubleFeld(A, Mphi, Mphi);
//	deleteDoubleFeld(b, Mphi);
//
//	return x;
//}
//
////
double* LGS_mit_alglib_loesen(double** A, double* b, int Mphi) {
//	real_1d_array bb;
//	real_2d_array AA;
//	AA.setlength(Mphi, Mphi);
//	bb.setlength(Mphi);
//
//	for (int i = 0; i < Mphi; ++i) {
//		bb[i] = b[i];
//		for (int k = 0; k < Mphi; ++k)
//			AA[i][k] = A[i][k];
//	}
//
//	real_1d_array e;
//	densesolverreport rep;
//	ae_int_t info;
//	rmatrixsolve(AA, Mphi, bb, info, rep, e);
//
//	double* erg = (double*) malloc(sizeof(double) * Mphi);
//	for (int k = 0; k < Mphi; ++k)
//		erg[k] = e[k];
//	return erg;
//	printf("Error 926\n");
	return NULL;
}

//
//int* bary(int k, int b) {
//     int* a;
//     if (k > 0) {
//         int jmax = (int) floor(log(k) / log(b));
//         a =(int*)malloc(sizeof(int)*(jmax + 1));
//         int q = (int) pow(b, jmax);
//         for (int j = 1; j <= jmax + 1; ++j) {
//             a[j - 1] = (int) floor((double) k / (double) q);
//             k = k - q * a[j - 1];
//             q = q / b;
//         }
//     }
//     return a;
// }
//
// int* nextbary(int l, int* ain, int b) {
//     int m = l;
//     int* aout =(int*)malloc(sizeof(int)*m);
//     bool carry = true;
//     for (int i = m; i >= 1; --i) {
//         if (carry) {
//             if (ain[i - 1] == b - 1) {
//                 aout[i - 1] = 0;
//             } else {
//                 aout[i - 1] = ain[i - 1] + 1;
//                 carry = false;
//             }
//         } else {
//             aout[i - 1] = ain[i - 1];
//         }
//     }
//     if (carry) {
//         int* Aout =(int*)malloc(sizeof(int)*(m + 1));
//         for (int i = 1; i < m+1; ++i) {
//             Aout[i] = aout[i - 1];
//         }
//         Aout[0] = 1;
//         return Aout;
//     } else {
//         return aout;
//     }
// }
//
// double** fauremat(int r, int i) {
//        double** C = DoubleFeld(r,r);
//        C[0][0] = 1;
//        for (int m = 2; m <= r; ++m) {
//            C[m - 1][m - 1] = 1;
//            C[0][m - 1] = i * C[0][m - 2];
//        }
//        for (int n = 3; n <= r; ++n) {
//            for (int m = 2; m < n; ++m) {
//                C[m - 1][n - 1] = C[m - 2][n - 2] + i * C[m - 1][n - 2];
//            }
//        }
//        return C;
//    }
//
// double mod(double a, double b)
//  {
//  int result = static_cast<int>( a / b );
//  return a - static_cast<double>( result ) * b;
//  }
//
//double** faurepts(int n0, int npts, int d, int b) {
//        int nmax = n0 + npts - 1;
//        int rmax = (int) (1 + floor(log(nmax) / log(b)));
//        double** P = DoubleFeld(npts,d);
//        for (int i = 0; i < npts; ++i) {
//            for (int k = 0; k < d; ++k) {
//                P[i][k] = 0;
//            }
//        }
//        double* y = DoubleFeld(rmax);
//        for (int i = 0; i < rmax; ++i) {
//            y[i] = 0;
//        }
//        int r = (int) (1 + floor(log(max(1, n0 - 1)) / log(b)));
//        double qnext = pow(b, r);
//        int* a = bary(n0 - 1, b);
//        int jmax = (int) floor(log(n0-1) / log(b));
//
//        double *** C = (double***)malloc(sizeof(double**)*d);
//
//        for (int i = 1; i < d; ++i) {
//            C[i] = fauremat(rmax, i);
//        }
//
//        double* bpwrs = DoubleFeld(rmax);
//        for (int i = 0; i < rmax; ++i) {
//            bpwrs[i] = 1. / pow(b, i + 1);
//        }
//        for (int k = n0; k <= nmax; ++k) {
//            a = nextbary(jmax+1,a, b);
//            if (k == qnext) {
//                r++;
//                qnext = b * qnext;
//            }
//            for (int j = 1; j <= r; ++j) {
//                P[k - n0][0] = P[k - n0][0] + bpwrs[j - 1] * a[r - j];
//            }
//            for (int i = 2; i <= d; ++i) {
//                for (int m = 1; m <= r; m++) {
//                    for (int n = 1; n <= r; ++n) {
//                        y[m - 1] = y[m - 1] + C[i - 1][m - 1][n - 1] * a[r - n];
//                    }
//                   // printf("vorher: %f\n",y[m-1]);
//                    y[m - 1] = y[m - 1]  -floor(y[m-1]/(double)b)*b;
//                 //   y[m - 1] = mod( y[m - 1]  ,b);
//                 //   printf("nachher: %f\n",y[m-1]);
//                    P[k - n0][i - 1] = P[k - n0][i - 1] + bpwrs[m - 1] * y[m - 1];
//                    y[m - 1] = 0;
//                }
//            }
//        }
//
//        for(int n=0;n<npts;n++)
//        	for(int k=0;k<d;++k)
//
//        		P[n][k]= P[n][k]  -floor(P[n][k]);
//        return P;
//    }

//bool schonfertig;

void tausche(double* daten, int * reihe, int i, int k) {
	double temp = daten[i];
	daten[i] = daten[k];
	daten[k] = temp;

	int temp2 = reihe[i];
	reihe[i] = reihe[k];
	reihe[k] = temp2;
}

int teile(double* daten, int* reihe, int links, int rechts) {
	int i = links;
	int j = rechts - 1;
	double pivot = daten[rechts];
	//if(pivot==10000)schonfertig=true;
	do {
		while ((daten[i] <= pivot) && i < rechts)
			i = i + 1;
		while ((pivot <= daten[j]) && j > links)
			j = j - 1;
		if (i < j)
			tausche(daten, reihe, i, j);
	} while (i < j);
	if (!(daten[i] <= pivot))
		tausche(daten, reihe, i, rechts);
	return i;
}

void quicksort(double* daten, int * reihe, int links, int rechts) {
	if (links < rechts) {
		int teiler = teile(daten, reihe, links, rechts);
		//	if(!schonfertig)
		quicksort(daten, reihe, links, teiler - 1);
		quicksort(daten, reihe, teiler + 1, rechts);
	}
}

int* quicksort(double* daten, int* reihe, int l) {
	int* reihe2 = new int[l];
	for (int ll = 0; ll < l; ++ll)
		reihe2[ll] = reihe[ll];

	double* daten2 = new double[l];

	for (int ll = 0; ll < l; ++ll)
		daten2[ll] = daten[reihe[ll]];

	quicksort(daten2, reihe2, 0, l - 1);
	delete[] daten2;

	//	int* reihe3=new int[l];
	//	for(int ll=0;ll<l;++ll)
	//			reihe3[ll]=reihe2[ll];
	//
	//delete[] reihe2;

	return reihe2;
}

int* quicksortK(double* daten, int l) {
	int* reihe = new int[l];
	for (int ll = 0; ll < l; ++ll)
		reihe[ll] = ll;

	quicksort(daten, reihe, 0, l - 1);
	return reihe;
}

int* quicksortStochPart(double* daten, int l, int number) {
	double wert = 0;
	double wertq = 0;
	for (int k = 0; k < 100; ++k) {
		wert += daten[k * (l - 1) / 100] / 100.;
		wertq += pow(daten[k * (l - 1) / 100] / 100., 2);
	}
	double schranke = wert - (0.6) * sqrt(pow(wert, 2) - wertq);

	int k = 0;
	for (int ll = 0; ll < l; ++ll)
		if (daten[ll] < schranke)
			k++;

	if (!(k > number)) {
		//		printf("Lieber normal \n");
		int* r = quicksort(daten, l);
		int* erg = new int[number];
		for (int n = 0; n < number; ++n)
			erg[n] = r[n];
		delete[] r;
		return erg;
	}

	double pre[k];
	int pos[k];
	k = 0;
	for (int ll = 0; ll < l; ++ll)
		if (daten[ll] < schranke) {
			pre[k] = daten[ll];
			pos[k] = ll;
			k++;
		}
	//int pin=rand();
	//printf("\n");
	//	for(int kk=0;kk<k;++kk)
	//	{		printf("%d:pos [%d]= %d\n",pin,kk,pos[kk]);
	//
	//		printf("%d:pos [%d]= %f\n",pin,kk,pre[kk]);
	//}
	int* r = quicksort((double*) pre, k);

	int* reihe = new int[number];
	for (int ll = 0; ll < number; ++ll)
		reihe[ll] = pos[r[ll]];

	delete[] r;

	return reihe;
}

int* quicksort(double* daten, int l) {
	int* reihe = new int[l];
	for (int ll = 0; ll < l; ++ll)
		reihe[ll] = ll;

	int* erg = quicksort(daten, reihe, l);
	delete[] reihe;
	return erg;
}

void quicksortUP(double* daten, int * reihe, int links, int rechts) {
	if (links < rechts) {
		int teiler = teile(daten, reihe, links, rechts);
		quicksort(daten, reihe, links, teiler - 1);
		if (abs(links - rechts) < 160)
			quicksort(daten, reihe, teiler + 1, rechts);
	}
}

int* quicksortUP(double* daten, int* reihe, int l) {
	int* reihe2 = new int[l];
	for (int ll = 0; ll < l; ++ll)
		reihe2[ll] = reihe[ll];

	double* daten2 = new double[l];

	for (int ll = 0; ll < l; ++ll)
		daten2[ll] = daten[reihe[ll]];

	quicksortUP(daten2, reihe2, 0, l - 1);
	delete[] daten2;

	//	int* reihe3=new int[l];
	//	for(int ll=0;ll<l;++ll)
	//			reihe3[ll]=reihe2[ll];
	//
	//delete[] reihe2;

	return reihe2;
}

int* quicksortUP(double* daten, int l) {
	int* reihe = new int[l];
	for (int ll = 0; ll < l; ++ll)
		reihe[ll] = ll;

	int* erg = quicksortUP(daten, reihe, l);
	delete[] reihe;
	return erg;
}

bool compare(const Such & a, const Such & b) {
	return a.wert < b.wert;
}

int* sort(vector<double> vec) {
	return partialsort(vec, vec.size());
}

int* partialsort(vector<double> vec, int n) {
	vector<Such> v;

	for (unsigned int l = 0; l < vec.size(); l++) {
		Such* s = new Such();
		s->pos = l;
		s->wert = vec.at(l);
		v.push_back(*s);
	}

	partial_sort(v.begin(), v.end(), v.begin() + n, compare);

	int* reihe = new int[n];
	for (int l = 0; l < n; ++l) {
		//		printf("%d, %f\n",v.at(l).pos,v.at(l).wert);
		reihe[l] = v.at(l).pos;
	}
	return reihe;
}
