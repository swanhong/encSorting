#ifndef PLAINSORTING_H_
#define PLAINSORTING_H_

#include "CyclicArray.h"
#include "vector"
#include "iostream"
#include "../MaskingGenerator.h"
#include "../PrintUtils.h"


class PlainSort {
public:
    MaskingGenerator* mg;
    long log2n;

    bool maskGen = false;
    bool maskTableGen = false;
    double*** mask;
    double*** maskRight;
    double*** maskTable;
    double*** maskTableRight;

    long logDataNum = 0;
    long colNum = 0;
    bool increase;
    long inc;


    PlainSort() {}
    PlainSort(long log2n);
    PlainSort(long log2n, long _logDataNum, long _colNum);
    ~PlainSort() {}

    void genMask();
    void genTableMask();
    void setInc(bool increase);

    void runPlainSorting(CyclicArray& ca, bool = true);

    long sortingRec(CyclicArray& ca, long logNum, long logJump, long loc);

    void compAndSwap(CyclicArray& ca, long loc, long dist);
    void compAndSwap(CyclicArray& ca, long loc, long dist, bool increase);

    void selfBitonicMerge(CyclicArray& ca, long log2n, double** mask, bool increase);

    void bitonicMerge(CyclicArray* ca, long log2n, long logNum);

    void bitonicMergeRec(CyclicArray* ca, long log2n, long start, long logNum, double** mask, double** mask2, bool increase);

    void compAndSwapTable(CyclicArray& ca, long loc, long dist, bool increase);
    void compAndSwapTable(CyclicArray& ca, long loc, long dist);
    
    void runPlainTableSorting(CyclicArray& ca, bool = true);

    long plainSortingTableRecursion(CyclicArray& ca, long logNum, long logJump, long loc);

    void bitonicTableMerge(CyclicArray* ca, long log2n, long logNum, long logDataNum, long colNum);
    
    void bitonicTableMergeRec(CyclicArray* ca, long log2n, long start, long logNum, long logDataNum, long colNum, double** maskCol, double*** maskInc, double*** maskDec, bool increase);

    void minMaxTable(CyclicArray& caLeft, CyclicArray& caRight, CyclicArray& caTableLeft, CyclicArray& caTableRight, long logDataNum, long colNum, double* maskRight, double* maskTableRight);

    void selfBitonicTableMerge(CyclicArray& ca, long log2n, long logDataNum, long colNum, double*** mask, bool increase);

};


#endif //! PLAINSORTING_H_
