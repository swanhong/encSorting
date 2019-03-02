#ifndef SORTINGALROTIGHM_H_
#define SORTINGALROTIGHM_H_

#include "CyclicArray.h"
#include "vector"
#include "iostream"
using namespace std;


class SortingAlgorithm {
public:
    CyclicArray ca;
    long log2n;
    long maskNum;
    double** mask;

    SortingAlgorithm(CyclicArray ca_, long log2n_);

    // void CompAndSwap(long start, long dist, long next, long num);
    void compAndSwap(long loc, long dist);
    void BatcherOddEvenSort();    
    long BatcherOddEvenSortRec(long logNum, long logJump, long loc);    
};

// void CompAndSwapOddEven(CyclicArray ca, long start, long dist, long next, long num, long log2n);

double* genMaskingComp(long length, long jump);
double* genMaskingMerge(long length, long num, long jump);

/*
    Generate all masking for Batch OESort
    length of mask = (log2n + 1) * log2n / 2
*/
void genAllMasking(long log2n, double**& mask);
long genMaskingRec(long log2n, long logNum, long logJump, double**& mask, long loc);



#endif // !SORTINGALROTIGHM_H_
