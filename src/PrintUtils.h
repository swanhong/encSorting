#ifndef PRINTUTILS_H_
#define PRINTUTILS_H_

#include "../HEAAN/src/HEAAN.h"
#include "Parameter.h"
#include "iostream"
#include <complex>
#include <math.h>
#include <string>

static bool WANT_TO_PRINT = true;
static bool DEC_AND_PRINT = false;

class PrintUtils {
public:
    static void parameter(Parameter param, std::string str);
    
    static void averageDifference(double* a1, std::complex<double>* a2, long n);
    static void averageDifference(std::complex<double>* a1, std::complex<double>* a2, long n);

    static void printArrays(double* a1, std::complex<double>* a2, long n);
    static void printArrays(std::complex<double>* a1, std::complex<double>* a2, long n);

    static void decAndPrint(std::string str, Ciphertext& cipher, Scheme& scheme, SecretKey& secretKey);

    static void decAndPrintTwo(std::string str, Ciphertext& cipher1, Ciphertext& cipher2, Scheme& scheme, SecretKey& secretKey);

    static void nprint(std::string str, bool isPrint);
};

;



#endif // !PRINTUTILS_H_