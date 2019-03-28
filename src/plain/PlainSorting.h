#ifndef PLAINSORTING_H_
#define PLAINSORTING_H_

#include "CyclicArray.h"
#include "vector"
#include "iostream"
#include "../MaskingGenerator.h"


class PlainSort {
public:
    PlainSort() {}
    ~PlainSort() {}
    void runPlainSorting(CyclicArray& ca, long log2n, bool=true);

    long sortingRec(CyclicArray& ca, double** mask, long logNum, long logJump, long loc, bool increase);

    void compAndSwap(CyclicArray& ca, double** mask, long loc, long dist, bool increase);
};


#endif //! PLAINSORTING_H_
