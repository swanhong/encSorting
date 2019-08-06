#ifndef ENCINTERFACE_H_
#define ENCINTERFACE_H_

#include "../../HEAAN/src/HEAAN.h"
#include "BootScheme.h"
#include "../Parameter.h"
#include "../PrintUtils.h"

class EncInterface {
public:
    Parameter param;
    Ring* ring;
    SecretKey* secretKey;
    BootScheme* scheme;
    BootHelper* bootHelper;
    bool printCondition;

    EncInterface(Parameter _param, bool = false);

    Ciphertext encrypt(double* mvec);
    complex<double>* decrypt(Ciphertext& cipher);
    ZZ* encode(double* mask);
    ZZ* flipPoly(ZZ* poly);
    void bootstrapping(Ciphertext& cipher);

    void add(Ciphertext& output, Ciphertext& a, Ciphertext& b);
    void sub(Ciphertext& output, Ciphertext& a, Ciphertext& b);
    void mult(Ciphertext& output, Ciphertext& a, Ciphertext& b);

    void nprint(string str);
    void nprint(string str, Ciphertext& cipher);

    void showTotalCount();
    void showCurrentCountAndReset();
};

#endif // !ENCINTERFACE_H_