/*
 * Hilfsmittel.h
 *
 *  Created on: Feb 6, 2012
 *      Author: dickmann
 */

#include <stdlib.h>

#include <stdlib.h>
#include <iostream>
#include <vector>
#include <stdio.h>
#include <math.h>
#include <fstream>
#include <cstring>

#ifndef HILFSMITTEL_H_
#define HILFSMITTEL_H_
using namespace std;

int argMin(vector<double> v);
int argZweiter(vector<double> v);
int argDritter(vector<double> v);
int argMax(vector<double> v);

double Min(vector<double> v);
double Max(vector<double> v);

void MatrixAusgeben(double**  a, int D);
vector<vector<double>> MatrixMultiplizieren(vector<vector<double>> a,vector<vector<double>>  b,int D);
double betrag(double x);

void BubbleSort(double* werte, int* index, int l);

double varianz(vector<double> vec);
double mittelwert(vector<double> vec);

int binary(int number, int digit);

void ErgebnisAnhaengen(vector<double> d, char* filename);
int* sort(vector<double> vec);
int* partialsort(vector<double> vec, int n);

int* array_machen(int z);
//int hargMax(double* v, int l, int ll);
//int hMax(double* v, int l, int ll);
double* LGSloesen(double** A, double* b, int Mphi);
double* LGS_mit_alglib_loesen(double** A, double* b, int Mphi);


void ErgebnisAnhaengen(double d);
void ErgebnisAnhaengenML(double d);
double* alphasLaden(int K);
void alphasSchreiben(double* alpha,int K);
void werteSchreiben(double* w,int K, int N);
void ErgebnisAnhaengen(double d, char* filename);

//double qnorm(double p);

vector<int> IntFeld(int m);
vector<vector<int>> IntFeld(int m,int n);
vector<vector<vector<int>>> IntFeld(int m,int n,int o);
vector<double>  DoubleFeld(int m);
vector<vector<double>> DoubleFeld(int m,int n);
vector<vector<vector<double>>> DoubleFeld(int m, int n, int o);
vector<vector<vector<vector<double>>>> DoubleFeld(int m, int n, int o, int p);
vector<vector<vector<vector<vector<double>>>>> DoubleFeld(int m, int n, int o, int p, int q);

double** faurepts(int n0, int npts, int d, int b);
void ausgeben(double* x, int j);

void InPipeSchreiben(int* pipe, double wert );
double AusPipeLesen(int* pipe);

double max(double x, double y);

double CumulativeNormalDistribution(double x);

void tausche(double* daten,int * reihe, int i, int k);
int teile(double* daten,int* reihe,int links, int rechts);
void quicksort(double* daten,int * reihe, int links, int rechts);
int* quicksortStochPart(double* daten, int l, int number);
int* quicksort(double* daten, int l);
int* quicksort(double* daten, int* reihe,  int l);
void quicksortUP(double* daten,int * reihe, int links, int rechts);
int* quicksortUP(double* daten, int l);
int* quicksortUP(double* daten, int* reihe,  int l);

int* quicksortK(double* daten, int l);

struct Such{
	double wert;
	int pos;
};

int* sort(vector<double> vec);

#endif /* HILFSMITTEL_H_ */
