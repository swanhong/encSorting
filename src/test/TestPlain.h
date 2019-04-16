#ifndef TESTPLAIN_H_
#define TESTPLAIN_H_

#include "../plain/CyclicArray.h"
#include "../plain/PlainSorting.h"
#include "../HEAAN/src/EvaluatorUtils.h"

class TestPlain {
public:
    static void plainSort(long logn, bool=true);

    static void showMasking(long log2n, bool increase);

    static void showBitonicMergeMasking(long log2n);

    static void bitonicMerge(long log2n, long logNum);
};


#endif // !TESTPLAIN_H_