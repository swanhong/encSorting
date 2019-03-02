#ifndef TESTALGORITHM_H_
#define TESTALGORITHM_H_

#include "CyclicArray.h"
#include "SortingAlrogithm.h"
#include "../HEAAN/src/HEAAN.h"

// #include "../ImprovedHEAAN-master/HEAAN/app/improved_bootstrapping/new_bootstrapping.h"
// #include "../ImprovedHEAAN-master/HEAAN/src/HEAAN.h"

struct Parameter {
    long logN; long logQ;
    long logp; long logc;
    long log2n; long radix;
    long logq;
    long logT;
};

class TestAlgorithm {
public:
    static void testBootstrapping(Parameter parameter);

    static void testOpreationsAndBoot(Parameter parameter);

    static void testMult(long logq, long logp, long logn);

    // static void testCompAndSwap(long length);

    static void testPlainSort(long logn);

    static void testMasking(long log2n);

    static void testSqrt(Parameter parameter);

    static void testMaxMin(Parameter parameter);

    static void testEncCompAndSwap(Parameter parameter);

    static void testEncSorting(Parameter parameter);
};

void fcnDecryptAndPrint(string str, Ciphertext cipher, Scheme scheme, SecretKey secretKey);

void fcnEncComputeSqrt(Ciphertext& outCipher, const Ciphertext& inCipher, long logp, long d, Scheme scheme);

void fcnEncMaxMin(Ciphertext& maxCipher, Ciphertext& minCipher, Ciphertext& input1, Ciphertext& input2, long logp, long iter, Scheme scheme);

void fcnGenEncAllMasking(long log2n, double**& mask);

long fcnGenEncMaskingRec(long log2n, long logNum, long logJump, double**& mask, long loc);

double* fcnGenEncMaskingComp(long length, long jump);

double* fcnGenEncMaskingMerge(long length, long num, long jump);

void fcnEncSort(Ciphertext& sortedCipher, const Ciphertext& inCipher, long logp, long iter, Scheme scheme);

void fcnEncCompAndSwap(Ciphertext& outCipher, const Ciphertext& inCipher, double* mask, long dist, long logp, long iter, Scheme scheme);

// int fcnEncBatcherOddEvenSortRec(Ciphertext& cipher, long num, long logJump, long loc, double** mask, long iter, Scheme scheme, BootHelper boothelper, Parameter parameter);

// void fcnEncBatcherOddEvenSort(Ciphertext& sortedCipher, const Ciphertext& inCipher, double** mask, long iter, Scheme scheme, BootHelper boothelper, Parameter parameter);

#endif // !TESTALGORITHM_H_