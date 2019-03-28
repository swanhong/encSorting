#include "PlainSorting.h"

void PlainSort::runPlainSorting(CyclicArray& ca, long log2n, bool increase) {
    MaskingGenerator mg(log2n, increase);
    double** mask = mg.getMasking();
    sortingRec(ca, mask, log2n, 0, 0, increase);
}

long PlainSort::sortingRec(CyclicArray& ca, double** mask, long logNum, long logJump, long loc, bool increase) {
    if (logNum == 1) {
        compAndSwap(ca, mask, loc, 1 << logJump, increase);
        return loc + 1;
    } else {
        if (logJump == 0) {
            loc = sortingRec(ca, mask, logNum - 1, logJump, loc, increase);
        }
        loc = sortingRec(ca, mask, logNum - 1, logJump + 1, loc, increase);
        compAndSwap(ca, mask, loc, 1 << logJump, increase);
        return loc + 1;
    }
}

void PlainSort::compAndSwap(CyclicArray& ca, double** mask, long loc, long dist, bool increase) {
    long length = ca.length;
    CyclicArray maskCA(mask[loc], length);     
    
    CyclicArray dummy(length);
    mult(dummy, ca, maskCA);
    ca.sub(dummy);
    if (increase) {
        dummy.rightRotate(dist);
        getMinMax(dummy, ca);
        dummy.leftRotate(dist);
    } else {
        dummy.leftRotate(dist);
        getMinMax(dummy, ca);
        dummy.rightRotate(dist);
    }  
    ca.add(dummy);
}