#include "RNG.h"
#include "math.h"
#include <stdlib.h>
#include <stdio.h>
#include <time.h>

#include <sys/time.h>

RNG::RNG() {
    int r = rand();
	
    struct timeval tim;
    	gettimeofday(&tim, NULL);

    setSeed( r +(int)(tim.tv_usec)+ time(NULL));
}

RNG::RNG(const RNG& orig) {

}

RNG::~RNG() {
}

void RNG::setSeed(int seed) {
    mt.init_genrand(seed);
    // Ein Mersennetwister braucht bei schlechten seeds bis zu
    // 600 Durchlaeufe um "auf Temperatur" zu kommen.
    for (int lauf = 0; lauf < 600; ++lauf)
        nextGaussian();
}

double RNG::BoxMuller(double U1, double U2) {
    if (U1 == 0)U1 = 0.000001;
    double R = -2 * log(U1);
    double V = 2. * 3.1415926535 * U2;
    return sqrt(R) * cos(V);
}

double RNG::nextUnif() {
    return mt.random();
}

double RNG::nextGaussian() {
    return BoxMuller(nextUnif(), nextUnif());
}
