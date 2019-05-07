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
    // cout << "compAndSwap with dist = " << dist << ", inc = " << increase << endl;
    // for (int i = 0; i < ca.length; i++) {
    //     cout << mask[loc][i] << " ";
    // }cout << endl;
    
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
    double** maskIncrease = mg.getBitonicMergeMasking();
    MaskingGenerator mg2(log2n, false);
    double** maskDecrease = mg2.getBitonicMergeMasking();

    bitonicMergeRec(ca, log2n, 0, logNum, maskIncrease, maskDecrease, true);
}

void PlainSort::bitonicMergeRec(CyclicArray* ca, long log2n, long start, long logNum, double** maskIncrease, double** maskDecrease, bool increase) {
    if (logNum == 0) return;        
    
    // cout << "logNum = " << logNum << endl;

    bitonicMergeRec(ca, log2n, start, logNum - 1, maskIncrease, maskDecrease, true);
    bitonicMergeRec(ca, log2n, start + (1 << (logNum - 1)), logNum - 1, maskIncrease, maskDecrease, false);
    
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
                // cout << "minMax (" << start + left << ", " << start + right << "), " << increase << endl;
                getMinMax(ca[start+left], ca[start + right]);
            }
        }
    }
    
    for(int i = 0; i < (1 << logNum); i++) {
        // cout << "self " << start + i << ", " << increase << endl;
        double ** mask;
        if (increase) mask = maskIncrease;
        else          mask = maskDecrease;
        
        selfBitonicMerge(ca[start + i], log2n, mask, increase);
    } 
    // cout << " -- end " << logNum << " -- " << endl;
}

void PlainSort::runPlainTableSorting(CyclicArray& ca, long log2n, long logDataNum, long colNum, bool increase) {
    MaskingGenerator mg(log2n, logDataNum, increase);
    double** mask = mg.getMasking();
    MaskingGenerator mg2(log2n, logDataNum, increase);
    double** maskOther = mg2.getMaskingOther();
    MaskingGenerator mgTable(log2n, logDataNum, colNum, increase);
    double** maskTable = mgTable.getMasking();
    MaskingGenerator mgTable2(log2n, logDataNum, colNum, increase);
    double** maskTableOther = mgTable2.getMaskingOther();

    long logn = log2n - logDataNum;
    long n = 1 << log2n;
    long maskNum = logn * (logn + 1) / 2;
    PrintUtils::printSingleMatrix("mask", mask, maskNum, n);
    PrintUtils::printSingleMatrix("maskOther", maskOther, maskNum, n);
    PrintUtils::printSingleMatrix("maskTable", maskTable, maskNum, n);
    PrintUtils::printSingleMatrix("maskTableOther", maskTableOther, maskNum, n);



    plainSortingTableRecursion(ca, logDataNum, colNum, log2n - logDataNum, 0, 0, mask, maskOther, maskTable, maskTableOther, increase);
}

long PlainSort::plainSortingTableRecursion(CyclicArray& ca, long logDataNum, long colNum, long logNum, long logJump, long loc, double** mask, double** maskOther, double** maskTable, double** maskTableOther, bool increase) {
    if (logNum == 1) {
        compAndSwapTable(ca, logDataNum, colNum, mask[loc], maskOther[loc], maskTable[loc], maskTableOther[loc], 1 << (logJump + logDataNum), increase);
    } else {
        if (logJump == 0) {
            loc = plainSortingTableRecursion(ca, logDataNum, colNum, logNum - 1, logJump, loc, mask, maskOther, maskTable, maskTableOther, increase);
        }
        loc = plainSortingTableRecursion(ca, logDataNum, colNum, logNum - 1, logJump + 1, loc, mask, maskOther, maskTable, maskTableOther, increase);
        compAndSwapTable(ca, logDataNum, colNum, mask[loc], maskOther[loc], maskTable[loc], maskTableOther[loc], 1 << (logJump + logDataNum), increase);
    }
    return loc + 1;
}

void PlainSort::compAndSwapTable(CyclicArray& ca, long logDataNum, long colNum, double* mask, double* maskRight, double* maskTable, double* maskTableRight, long dist, bool increase) {
    long length = ca.length;
    CyclicArray maskPoly(mask, length);
    CyclicArray maskRightPoly(maskRight, length);
    CyclicArray maskTablePoly(maskTable, length);
    CyclicArray maskTableRightPoly(maskTableRight, length);
    CyclicArray ca1, ca1Right, caTable, caTableRight;
    mult(ca1, ca, maskPoly);
    mult(ca1Right, ca, maskRightPoly);
    mult(caTable, ca, maskTablePoly);
    mult(caTableRight, ca, maskTableRightPoly);

    // ca.printAsVector();
    // maskTablePoly.printAsVector();
    // caTable.printAsVector();

    ca.sub(ca1);
    ca.sub(ca1Right);

    if(increase) caTable.rightRotate(dist);
    else caTable.leftRotate(dist);

    cout << "before comp" << endl;
    caTable.printAsVector();
    caTableRight.printAsVector();

    for (int i = 0; i < length; i++) {
        if(caTable.get(i) >= caTableRight.get(i)) {
            caTable.set(i, 1);
            caTableRight.set(i, 0);
        } else {
            caTable.set(i, 0);
            caTableRight.set(i, 1);
        }
    }

    cout << "after comp" << endl;
    caTable.printAsVector();
    caTableRight.printAsVector();
    maskTableRightPoly.printAsVector();
    
    caTable.mult(maskTableRightPoly);
    caTableRight.mult(maskTableRightPoly);
    
    if(increase) caTable.leftRotate(dist);
    else caTable.rightRotate(dist);

    cout << "mult mask and rot" << endl;
    caTable.printAsVector();
    caTableRight.printAsVector();

    caTable.leftRotate(colNum);
    caTableRight.leftRotate(colNum);

    for (int i = 0; i < logDataNum; i++) {
        CyclicArray tmpTable(caTable);
        CyclicArray tmpTableRight(caTableRight);
    
        tmpTable.rightRotate(1 << i);
        tmpTableRight.rightRotate(1 << i);
    
        caTable.add(tmpTable);
        caTableRight.add(tmpTableRight);
    }

    cout << "copy" << endl;
    caTable.printAsVector();
    caTableRight.printAsVector();
    
    CyclicArray caTableFlip(caTable);
    CyclicArray caTableFlipRight(caTableRight);
    for (int i = 0; i < length; i++) {
        caTableFlip.data[i] *= -1;
        caTableFlipRight.data[i] *= -1;
    }
    caTableFlip.add(maskPoly);
    caTableFlipRight.add(maskRightPoly);

    cout << "flip" << endl;
    caTableFlip.printAsVector();
    caTableFlipRight.printAsVector();
    cout << endl;

    CyclicArray caLeftSmall, caLeftBig, caRightSmall, caRightBig;
    mult(caLeftSmall, ca1, caTableFlip);
    mult(caLeftBig, ca1, caTable);
    mult(caRightSmall, ca1Right, caTableFlipRight);
    mult(caRightBig, ca1Right, caTableRight);

    caLeftSmall.printAsVector();
    caLeftBig.printAsVector();
    caRightSmall.printAsVector();
    caRightBig.printAsVector();

    if (increase) {
        caLeftBig.rightRotate(dist);
        caRightSmall.leftRotate(dist);
    } else {
        caLeftBig.leftRotate(dist);
        caRightSmall.rightRotate(dist);
    }
    
    ca.add(caLeftSmall);
    ca.add(caLeftBig);
    ca.add(caRightSmall);
    ca.add(caRightBig);

    // cout << " ===========================" << endl;
    // double* mvec = ca.getArray();
    // PrintUtils::printArraysWithDataNum(mvec, mvec, length, logDataNum, 0);
}