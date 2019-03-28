#ifndef TESTPLAIN_H_
#define TESTPLAIN_H_

#include "../plain/CyclicArray.h"
#include "../plain/PlainSorting.h"

class TestPlain {
public:
    static void plainSort(long logn, bool increase);

    static void showMasking(long log2n, bool increase);

    static void showBitonicMergeMasking(long log2n);
};


#endif // !TESTPLAIN_H_