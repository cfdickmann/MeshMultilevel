#ifndef LCG_GENERATOR_H
#define	LCG_GENERATOR_H

#include "mt.h"
#include <iostream>
#include <math.h>
#include <stdlib.h>
#include "math.h"
#include "stdlib.h"
#include "Hilfsmittel.h"

class RNG {
public:
    RNG();
    RNG(const RNG& orig);
    virtual ~RNG();
    MersenneTwister mt;
    void setSeed(int seed);
    double nextUnif();
	double BoxMuller(double U1, double U2);
    double nextGaussian();

    void generate_numbers();
};

#endif	/* LCG_GENERATOR_H */

