#ifndef TESTALGORITHM_H_
#define TESTALGORITHM_H_

#include "EncContext.h"
#include "CyclicArray.h"
#include "SortingAlgorithm.h"
#include "Parameter.h"

class TestAlgorithm {
public:

    static void testBootstrapping(Parameter parameter);

    static void testOpreationsAndBoot(Parameter parameter);

    static void testMult(Parameter parameter);

    // static void testCompAndSwap(long length);

    static void testPlainSort(long logn);

    static void testMasking(long log2n);

    static void testSqrt(Parameter parameter);

    static void testMaxMin(Parameter parameter);

    static void testEncCompAndSwap(Parameter parameter, long iter);
};

double difference(double* a1, complex<double>* a2, long n);

// int fcnEncBatcherOddEvenSortRec(Ciphertext& cipher, long num, long logJump, long loc, double** mask, long iter, Scheme scheme, BootHelper boothelper, Parameter parameter);

// void fcnEncBatcherOddEvenSort(Ciphertext& sortedCipher, const Ciphertext& inCipher, double** mask, long iter, Scheme scheme, BootHelper boothelper, Parameter parameter);

#endif // !TESTALGORITHM_H_