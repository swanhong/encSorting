#include "PlainSorting.h"

void PlainSort::runSorting(CyclicArray& ca, long log2n) {
    MaskingGenerator mg(log2n);
    double** mask = mg.getMasking();
    sortingRec(ca, mask, log2n, 0, 0);
}

long PlainSort::sortingRec(CyclicArray& ca, double** mask, long logNum, long logJump, long loc) {
    if (logNum == 1) {
        compAndSwap(ca, mask, loc, 1 << logJump);
        return loc + 1;
    } else {
        if (logJump == 0) {
            loc = sortingRec(ca, mask, logNum - 1, logJump, loc);
        }
        loc = sortingRec(ca, mask, logNum - 1, logJump + 1, loc);
        compAndSwap(ca, mask, loc, 1 << logJump);
        return loc + 1;
    }
}

void PlainSort::compAndSwap(CyclicArray& ca, double** mask, long loc, long dist) {
    long length = ca.length;
    CyclicArray maskCA(mask[loc], length);     
    
    CyclicArray dummy(length);
    mult(dummy, ca, maskCA);
    ca.sub(dummy);
    dummy.rightRotate(dist);
    getMinMax(dummy, ca);
    dummy.leftRotate(dist);
    ca.add(dummy);
}