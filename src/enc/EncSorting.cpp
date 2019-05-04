#include "EncSorting.h"


void EncSorting::runEncSorting(Ciphertext& cipher, BootScheme& scheme, Ring& ring, BootHelper& bootHelper, bool increase) {
        MaskingGenerator mg(param.log2n, increase);
        double** mask = mg.getMasking();
        bootAlgo = BootAlgo(param, sqrtIter, increase);
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

void EncSorting::runEncSortingDec(Ciphertext& cipher, BootScheme& scheme, Ring& ring, BootHelper& bootHelper, bool increase, SecretKey sk) {
        MaskingGenerator mg(param.log2n, increase);
        double** mask = mg.getMasking();
        bootAlgo = BootAlgo(param, sqrtIter, increase);
        PlainSort plainSort;
        sortingRecursion(cipher, param.log2n, 0, 0, mask, scheme, ring, bootHelper, sk, plainSort);
        scheme.showTotalCount();
}

long EncSorting::sortingRecursion(Ciphertext& cipher, long logNum, long logJump, long loc, double** mask, BootScheme& scheme, Ring& ring, BootHelper& bootHelper, SecretKey sk, PlainSort ps) {
    TimeUtils timeutils;
    PrintUtils::nprint("run encBatcherSort with loc = " + to_string(loc), WANT_TO_PRINT);
    if (logNum == 1) {
        timeutils.start(to_string(loc)+"th CompAndSwap");
        
        complex<double>* mvecCplx = scheme.decrypt(sk, cipher);
        double* mvec = new double[cipher.n];
        for (int i = 0; i < cipher.n; i++) mvec[i] = mvecCplx[i].real();
        CyclicArray ca(mvec, cipher.n);
        ps.compAndSwap(ca, mask, loc, 1 << logJump, true);
        mvec = ca.getArray();
        
        bootAlgo.compAndSwapDec(cipher, mask[loc], 1<< logJump, scheme, ring, bootHelper, loc, sk);
        
        complex<double>* dvec = scheme.decrypt(sk, cipher);
        cout << loc << "th compAndSwap Result" << endl;
        PrintUtils::printArrays(mvec, dvec, cipher.n);
        cout << endl << "******************" << endl;
        PrintUtils::averageDifference(mvec, dvec, cipher.n);
        cout << "******************" << endl << endl;

        scheme.showCurrentCount();
        scheme.resetCount();
        timeutils.stop(to_string(loc)+"th CompAndSwap");
    } else {
        if (logJump == 0) {
            loc = sortingRecursion(cipher, logNum - 1, logJump, loc, mask, scheme, ring, bootHelper, sk, ps);
        }
        loc = sortingRecursion(cipher, logNum - 1, logJump + 1, loc, mask, scheme, ring, bootHelper, sk, ps);

        timeutils.start(to_string(loc)+"th CompAndSwap");

        complex<double>* mvecCplx = scheme.decrypt(sk, cipher);
        double* mvec = new double[cipher.n];
        for (int i = 0; i < cipher.n; i++) mvec[i] = mvecCplx[i].real();
        CyclicArray ca(mvec, cipher.n);
        ps.compAndSwap(ca, mask, loc, 1 << logJump, true);
        mvec = ca.getArray();

        bootAlgo.compAndSwapDec(cipher, mask[loc], 1<< logJump, scheme, ring, bootHelper, loc, sk);

        complex<double>* dvec = scheme.decrypt(sk, cipher);
        cout << loc << "th compAndSwap Result" << endl;
        PrintUtils::printArrays(mvec, dvec, cipher.n);
        cout << endl << "******************" << endl;
        PrintUtils::averageDifference(mvec, dvec, cipher.n);
        cout << "******************" << endl << endl;

        scheme.showCurrentCount();
        scheme.resetCount();
        timeutils.stop(to_string(loc)+"th CompAndSwap");
    }
    return loc + 1;
}

void EncSorting::runEncTableSorting(Ciphertext& cipher, long logDataNum, long colNum, BootScheme& scheme, Ring& ring, BootHelper& bootHelper, SecretKey& secretKey, bool increase) {
    MaskingGenerator mg(param.log2n, logDataNum);
    double** mask = mg.getMasking();
    MaskingGenerator mg2(param.log2n, logDataNum);
    double** maskOther = mg2.getMaskingOther();
    MaskingGenerator mgTable(param.log2n, logDataNum, colNum, true);
    double** maskTable = mgTable.getMasking();
    MaskingGenerator mgTable2(param.log2n, logDataNum, colNum, true);
    double** maskTableOther = mgTable2.getMaskingOther();

    
    long logn = param.log2n - logDataNum;
    long n = 1 << param.log2n;
    long maskNum = logn * (logn + 1) / 2;
    PrintUtils::printSingleMatrix("mask", mask, maskNum, n);
    cout << endl;
    PrintUtils::printSingleMatrix("maskOther", maskOther, maskNum, n);
    cout << endl;
    PrintUtils::printSingleMatrix("maskTable", maskTable, maskNum, n);
    cout << endl;
    PrintUtils::printSingleMatrix("maskTableOther", maskTableOther, maskNum, n);

    bootAlgo = BootAlgo(param, invIter, compIter);
    PlainSort plainSort;
    sortingTableRecursion(cipher, logDataNum, param.log2n - logDataNum, 0, 0, mask, maskOther, maskTable, maskTableOther, scheme, ring, bootHelper, secretKey, plainSort);
    scheme.showTotalCount();
}

long EncSorting::sortingTableRecursion(Ciphertext& cipher, long logDataNum, long logNum, long logJump, long loc,
                                    double** mask, double** maskRight, double** maskTable, double** maskTableRight,
                                    BootScheme& scheme, Ring& ring, BootHelper& bootHelper, SecretKey& secretKey, PlainSort ps) {
    
    std::cout << "(" << logNum << ", " << logJump << ", " << loc << ")" << std::endl;
    TimeUtils timeutils;
    if (logNum == 1) {
        timeutils.start(to_string(loc)+"th CompAndSwap");

        CompAndSwapTableBothWithDec(cipher, logDataNum, 1 << (logJump + logDataNum), mask[loc], maskRight[loc], maskTable[loc], maskTableRight[loc], scheme, ring, bootHelper, secretKey, ps);

        scheme.showCurrentCount();
        scheme.resetCount();
        timeutils.stop(to_string(loc)+"th CompAndSwap with " + to_string(logNum) + ", " + to_string(logJump));
    } else {
        if (logJump == 0) {
            loc = sortingTableRecursion(cipher, logDataNum, logNum - 1, logJump, loc, mask, maskRight, maskTable, maskTableRight, scheme, ring, bootHelper, secretKey, ps);
        }
        loc = sortingTableRecursion(cipher, logDataNum, logNum - 1, logJump + 1, loc, mask, maskRight, maskTable, maskTableRight, scheme, ring, bootHelper, secretKey, ps);
        timeutils.start(to_string(loc)+"th CompAndSwap");

        CompAndSwapTableBothWithDec(cipher, logDataNum, 1 << (logJump + logDataNum), mask[loc], maskRight[loc], maskTable[loc], maskTableRight[loc], scheme, ring, bootHelper, secretKey, ps);
        scheme.showCurrentCount();
        scheme.resetCount();
        timeutils.stop(to_string(loc)+"th CompAndSwap with " + to_string(logNum) + ", " + to_string(logJump));
    }
    return loc + 1;
}

void EncSorting::CompAndSwapTableBothWithDec(Ciphertext& cipher, long logDataNum, long dist,
                                    double* mask, double* maskRight, double* maskTable, double* maskTableRight,
                                    BootScheme& scheme, Ring& ring, BootHelper& bootHelper, SecretKey& secretKey, PlainSort ps) {
        cout << "run CompAndSwapTableBothWithDec with logDataNum = " << logDataNum << ", dist = " << dist << endl;
        complex<double>* mvecCplx = scheme.decrypt(secretKey, cipher);
        double* mvec = new double[cipher.n];
        for (int i = 0; i < cipher.n; i++) mvec[i] = mvecCplx[i].real();
        CyclicArray ca(mvec, cipher.n);
        
        long logq = cipher.logq;
        long logp = cipher.logp;
        bootAlgo.compAndSwapTable(cipher, logDataNum, mask, maskRight, maskTable, maskTableRight, dist, scheme, ring, bootHelper, secretKey);
        ps.compAndSwapTable(ca, logDataNum, mask, maskRight, maskTable, maskTableRight, dist);

        mvec = ca.getArray();
        complex<double>* dvec = scheme.decrypt(secretKey, cipher);
        cout << " compAndSwap Result" << endl;
        cout << "cipher.logq, logp = " << logq << ", " << logp << " -> " << cipher.logq << ", " << cipher.logp << endl;
        PrintUtils::printArrays(mvec, dvec, cipher.n);
        cout << endl << "******************" << endl;
        PrintUtils::averageDifference(mvec, dvec, cipher.n);
        cout << "******************" << endl << endl;
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

    BootAlgo bootAlgo(param, sqrtIter, increase);

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