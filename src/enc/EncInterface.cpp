#include "EncInterface.h"

EncInterface::EncInterface(Parameter _param, bool _printCondition) :
    param(_param), printCondition(_printCondition) {
    srand(time(NULL));
    SetNumThreads(16);
    long n = 1 << log2n;
    
    PrintUtils::parameter(param, "HEAAN");

    ring = new Ring(logN, logQ);
    secretKey = new SecretKey(*ring);
    scheme = new BootScheme(*secretKey, *ring);
    scheme->addConjKey(*secretKey);
    scheme->addLeftRotKeys(*secretKey);
    scheme->addRightRotKeys(*secretKey);
    bootHelper = new BootHelper(log2n, radix, logc, *scheme, *ring, *secretKey);
}

Ciphertext EncInterface::encrypt(double* mvec) {
    return scheme->encrypt(mvec, 1 << log2n, logp, logQ);
}

complex<double>* EncInterface::decrypt(Ciphertext& cipher) {
    return scheme->decrypt(*secretKey, cipher);
}

ZZ* EncInterface::encode(double* mask) {
    ZZ* maskPoly = new ZZ[1 << logN];
    ring->encode(maskPoly, mask, 1 << log2n, logp);
    return maskPoly;
}

ZZ* EncInterface::flipPoly(ZZ* poly) {
    ZZ* res = new ZZ[1 << logN];
    for (int i = 0; i < (1 << logN); i++) {
        res[i] = -poly[i];
    }
    res[0] += 1;
    return res;    
}

void EncInterface::bootstrapping(Ciphertext& cipher) {
    bootHelper->bootstrapping_cos(cipher, logq, logQ, 5);
}

void EncInterface::add(Ciphertext& output, Ciphertext& a, Ciphertext& b) {
    output = scheme->add(a, b);
}

void EncInterface::sub(Ciphertext& output, Ciphertext& a, Ciphertext& b) {
    output = scheme->sub(a, b);
}

void EncInterface::mult(Ciphertext& output, Ciphertext& a, Ciphertext& b) {
    output = scheme->mult(a, b);
    scheme->reScaleByAndEqual(output, logp);
    scheme->resetImagErrorAndEqual(output);
}

void EncInterface::nprint(string str) {
    if(printCondition) {
        cout << str << endl;
    }
}

void EncInterface::nprint(string str, Ciphertext& cipher) {
    if(printCondition) {
        scheme->decryptAndPrint(str, *secretKey, cipher);
        std::cout << "logq = " << cipher.logq << endl;
    }
}

void EncInterface::showTotalCount() {
    scheme->showTotalCount();
}

void EncInterface::showCurrentCountAndReset() {
    scheme->showCurrentCount();
    scheme->resetCount();
}