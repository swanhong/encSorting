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

void PlainSort::selfBitonicMerge(CyclicArray& ca, long log2n, double** mask, bool increase) {
    for(int i = 0; i <log2n; i++) {
        compAndSwap(ca, mask, i, 1 << (log2n - 1 - i), increase);
    }
}

void PlainSort::bitonicMerge(CyclicArray* ca, long log2n, long logNum) {
    MaskingGenerator mg(log2n);
    double** mask = mg.getBitonicMergeMasking();
    MaskingGenerator mg2(log2n, false);
    double** mask2 = mg2.getBitonicMergeMasking();
    
    bitonicMergeRec(ca, log2n, 0, logNum, mask, mask2, true);
}

void PlainSort::bitonicMergeRec(CyclicArray* ca, long log2n, long start, long logNum, double** mask, double** mask2, bool increase) {
    if (logNum == 0) return;

    cout << "logNum = " << logNum << endl;

    bitonicMergeRec(ca, log2n, start, logNum - 1, mask, mask2, true);
    bitonicMergeRec(ca, log2n, start + (1 << (logNum - 1)), logNum - 1, mask, mask2, false);

    for(int i = 0; i < logNum; i++) {
        for(int j = 0; j < (1 << i); j++) {
            for(int k = 0; k < (1 << (logNum - 1 - i)); k++) {
                long left = 2 * j * (1 << (logNum - i - 1)) + k;
                long right = (2 * j + 1) * (1 << (logNum - i - 1)) + k;
                if (!increase) {
                    long x = left;
                    left = right;
                    right = x;
                }
                cout << "minMax ( " << start + left << ", " << start + right << "), " << increase << endl;
                getMinMax(ca[start+left], ca[start + right]);
            }
        }
    }
    for(int i = 0; i < (1 << logNum); i++) {
        cout << start + i << " : "; ca[start + i].printAsVector();
    }
    
    for(int i = 0; i < (1 << logNum); i++) {
        cout << "self " << start + i << ", " << increase << endl;
        if (increase) {
            selfBitonicMerge(ca[start + i], log2n, mask, true);
        } else {
            selfBitonicMerge(ca[start + i], log2n, mask2, false);
        }
    } 
    cout << "at logNum = " << logNum << endl;
    for(int i = 0; i < (1 << logNum); i++) {
        cout << start + i << " : "; ca[start + i].printAsVector();
    }
    
    cout << " -- end " << logNum << " -- " << endl;
}

void PlainSort::runPlainSorting(double* mvec, long log2n, bool increase) {
    long length = 1 << log2n;
    CyclicArray ca(mvec, length);
    runPlainSorting(ca, log2n, increase);
    ca.getArray(mvec);
}