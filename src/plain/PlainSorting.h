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
    void runSorting(CyclicArray& ca, long log2n);

    long sortingRec(CyclicArray& ca, double** mask, long logNum, long logJump, long loc);

    void compAndSwap(CyclicArray& ca, double** mask, long loc, long dist);
};


#endif //! PLAINSORTING_H_
