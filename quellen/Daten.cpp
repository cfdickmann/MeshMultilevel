#include "AmericanOption.h"

using namespace std;

void AmericanOption::Daten() {
    int Example =3;

    X0 = DoubleFeld(100); //genug Platz
    sigma = DoubleFeld(100);

//    if (Example == 3) { //Glasserman Example MaxCall
//        option = MIN_PUT;
//        delta = 0.0;
//        D = 2;//Achtung!
//        for (int j = 0; j < D; ++j) {
//            X0[j] = 100.;
//            sigma[j] = 0.2;
//        }
//        Strike = 100.;
//        r = 0.05;
//        T = 3;
//        N = 10;
//    }

    if (Example == 3) { //Glasserman Example MaxCall
            option = MAX_CALL;
            delta = 0.1;
            D = 2;//Achtung!
            for (int j = 0; j < D; ++j) {
                X0[j] = 9.;
                sigma[j] = 0.2;
            }
            Strike = 10.;
            r = 0.05;
            T = 3;
            N = 10;
        }

}