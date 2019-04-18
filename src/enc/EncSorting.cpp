#include "EncSorting.h"


void EncSorting::runEncSorting(Ciphertext& cipher, BootScheme& scheme, Ring& ring, BootHelper& bootHelper, bool increase) {
        MaskingGenerator mg(param.log2n, increase);
        double** mask = mg.getMasking();
        bootAlgo = BootAlgo(param, iter, increase);
        sortingRecursion(cipher, param.log2n, 0, 0, mask, scheme, ring, bootHelper);
        scheme.showTotalCount();
}

long EncSorting::sortingRecursion(Ciphertext& cipher, long logNum, long logJump, long loc, double** mask, BootScheme& scheme, Ring& ring, BootHelper& bootHelper) {
    TimeUtils timeutils;
    PrintUtils::nprint("run encBatcherSort with loc = " + to_string(loc), WANT_TO_PRINT);
    if (logNum == 1) {
        timeutils.start(to_string(loc)+"th CompAndSwap");
        bootAlgo.compAndSwap(cipher, mask[loc], 1<< logJump, scheme, ring, bootHelper);
        scheme.showCurrentCount();
        scheme.resetCount();
        timeutils.stop(to_string(loc)+"th CompAndSwap");
    } else {
        if (logJump == 0) {
            loc = sortingRecursion(cipher, logNum - 1, logJump, loc, mask, scheme, ring, bootHelper);
        }
        loc = sortingRecursion(cipher, logNum - 1, logJump + 1, loc, mask, scheme, ring, bootHelper);

        timeutils.start(to_string(loc)+"th CompAndSwap");
        bootAlgo.compAndSwap(cipher, mask[loc], 1<< logJump, scheme, ring, bootHelper);
        scheme.showCurrentCount();
        scheme.resetCount();
        timeutils.stop(to_string(loc)+"th CompAndSwap");
    }
    return loc + 1;
}

void EncSorting::bitonicMerge(Ciphertext* cipher, long logNum, BootScheme& scheme, Ring& ring, BootHelper& bootHelper) {
    MaskingGenerator mg(param.log2n);
    double** maskIncrease = mg.getBitonicMergeMasking();
    MaskingGenerator mg2(param.log2n, false);
    double** maskDecrease = mg2.getBitonicMergeMasking();

    bitonicMergeRec(cipher, 0, logNum, maskIncrease, maskDecrease, scheme, ring, bootHelper, true);

    scheme.showTotalCount();
}

void EncSorting::bitonicMergeRec(Ciphertext* cipher, long start, long logNum, double** maskIncrease, double** maskDecrease, BootScheme& scheme, Ring& ring, BootHelper& bootHelper, bool increase) {
    if (logNum == 0) return;        
    
    bitonicMergeRec(cipher, start, logNum - 1, maskIncrease, maskDecrease, scheme, ring, bootHelper, true);
    bitonicMergeRec(cipher, start + (1 << (logNum - 1)), logNum - 1, maskIncrease, maskDecrease, scheme, ring, bootHelper, false);

    BootAlgo bootAlgo(param, iter, increase);

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
                cout << "-- Run minMax (" << start + left << ", " << start + right << ")" << endl;
                bootAlgo.minMax(cipher[start + left], cipher[start + right], scheme ,bootHelper);
            }
        }
    }
    scheme.showCurrentCount();
    scheme.resetCount();
    
    double** mask;
    if (increase) mask = maskIncrease;
    else          mask = maskDecrease;

    cout << "----- run selfBitonicMerge " << endl;
    for(int i = 0; i < (1 << logNum); i++) {        
        cout << "           -- ctxt " << start + i << endl;
        bootAlgo.selfBitonicMerge(cipher[start + i], mask, scheme, ring, bootHelper);
    } 
    scheme.showCurrentCount();
    scheme.resetCount();
    cout << "-----";
}

void EncSorting::reverseHalf(Ciphertext* cipher, long logNum, BootScheme& scheme, Ring& ring, BootHelper& bootHelper) {
    long num = 1 << logNum;
    MaskingGenerator mg(param.log2n);
    double** mask = mg.getBitonicMergeMasking();    
    BootAlgo bootAlgo(param, 0);
    cout << "---- Reverse odd ctxt " << endl;
    for (int i = 0; i < num / 2; i++) {
        
        bootAlgo.reverse(cipher[2 * i + 1], mask, scheme, ring, bootHelper);
    }
    cout << "---- modDown even ctxt " << endl;
    for (int i = 0; i < num / 2; i++) {
        
        scheme.modDownToAndEqualModified(cipher[2 * i], cipher[2 * i + 1], bootHelper, param);
    }    
} 