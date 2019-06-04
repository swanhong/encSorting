#include "PlainSorting.h"


PlainSort::PlainSort(long _log2n) : log2n(_log2n) {
    genMask();
}

PlainSort::PlainSort(long _log2n, long _logDataNum, long _colNum) : 
    log2n(_log2n), logDataNum(_logDataNum), colNum(_colNum) {
    genTableMask();
}

void PlainSort::setInc(bool increase) {
    if (increase) {
        inc = 0;
    } else {
        inc = 1;
    }
}

void PlainSort::genMask() {
    if(!maskGen) {
        maskGen = true;
        mask = new double**[2];
        for (int inc = 0; inc < 2; inc++) {
            mg = new MaskingGenerator(log2n, (inc == 0));
            mask[inc] = mg->getMasking();
        }
    }
}

void PlainSort::genTableMask() {
    if(!maskTableGen) {
        maskTableGen = true;

        mask = new double**[2];
        for (int inc = 0; inc < 2; inc++) {
            mg = new MaskingGenerator(log2n, logDataNum, inc == 0);
            mask[inc] = mg->getMasking();
        }
        
        maskRight = new double**[2];
        for (int inc = 0; inc < 2; inc++) {
            mg = new MaskingGenerator(log2n, logDataNum, inc == 0);
            maskRight[inc] = mg->getMaskingOther();
        }

        maskTable = new double**[2];
        for (int inc = 0; inc < 2; inc++) {
            mg = new MaskingGenerator(log2n, logDataNum, colNum, inc == 0);
            maskTable[inc] = mg->getMasking();
        }

        maskTableRight = new double**[2];
        for (int inc = 0; inc < 2; inc++) {
            mg = new MaskingGenerator(log2n, logDataNum, colNum, inc == 0);
            maskTableRight[inc] = mg->getMaskingOther();
        }
        std::cout << "check mask" << std::endl;
        mg->printMask(maskTableRight[0], mg->maskNum);
    }
}

void PlainSort::genMergeMask() {
    if(!maskMergeGen) {
        maskMergeGen = true;
        maskMerge = new double**[2];
        for (int inc = 0; inc < 2; inc++) {
            mg = new MaskingGenerator(log2n, (inc == 0));
            mask[inc] = mg->getBitonicMergeMasking();
        }
    }
}

void PlainSort::runPlainSorting(CyclicArray& ca, bool _increase) {
    increase = _increase;
    setInc(_increase);
    genMask();
    sortingRec(ca, log2n, 0, 0);
}

long PlainSort::sortingRec(CyclicArray& ca, long logNum, long logDist, long loc) {
    if (logNum == 1) {
        compAndSwap(ca, loc, 1 << logDist);
    } else {
        if (logDist == 0) {
            loc = sortingRec(ca, logNum - 1, logDist, loc);
        }
        loc = sortingRec(ca, logNum - 1, logDist + 1, loc);
        compAndSwap(ca, loc, 1 << logDist);
    }
    return loc + 1;
}

void PlainSort::compAndSwap(CyclicArray& ca, long loc, long dist, bool _increase) {
    increase = _increase;
    setInc(increase);
    compAndSwap(ca, loc, dist);
}

void PlainSort::compAndSwap(CyclicArray& ca, long loc, long dist) {
    // cout << "compAndSwap with dist = " << dist << ", inc = " << increase << endl;
    // for (int i = 0; i < ca.length; i++) {
    //     cout << mask[loc][i] << " ";
    // }cout << endl;
    long length = ca.length;
    CyclicArray maskCA(mask[inc][loc], length);
    
    CyclicArray dummy(length);
    mult(dummy, ca, maskCA);

    ca.sub(dummy);
    
    dummy.rightRotateConditional(dist, increase);
    getMinMax(dummy, ca);
    dummy.leftRotateConditional(dist, increase);
    ca.add(dummy);
}

void PlainSort::selfBitonicMerge(CyclicArray& ca, long log2n, double** mask, bool increase) {
    for(int i = 0; i <log2n; i++) {
        compAndSwap(ca, i, 1 << (log2n - 1 - i), increase);
    }
}

void PlainSort::bitonicMerge(CyclicArray* ca, long log2n, long logNum) {
    genMergeMask();
    bitonicMergeRec(ca, log2n, 0, logNum, true);
}

void PlainSort::bitonicMergeRec(CyclicArray* ca, long log2n, long start, long logNum, bool increase) {
    if (logNum == 0) return;        
    
    // cout << "logNum = " << logNum << endl;

    bitonicMergeRec(ca, log2n, start, logNum - 1, true);
    bitonicMergeRec(ca, log2n, start + (1 << (logNum - 1)), logNum - 1, false);
    
    setInc(increase);

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
        selfBitonicMerge(ca[start + i], log2n, maskMerge[inc], increase);
    } 
    // cout << " -- end " << logNum << " -- " << endl;
}

void PlainSort::runPlainTableSorting(CyclicArray& ca, bool _increase) {
    increase = _increase;
    setInc(increase);
    // MaskingGenerator mg(log2n, logDataNum, increase);
    // double** mask = mg.getMasking();
    // MaskingGenerator mg2(log2n, logDataNum, increase);
    // double** maskOther = mg2.getMaskingOther();
    // MaskingGenerator mgTable(log2n, logDataNum, colNum, increase);
    // double** maskTable = mgTable.getMasking();
    // MaskingGenerator mgTable2(log2n, logDataNum, colNum, increase);
    // double** maskTableOther = mgTable2.getMaskingOther();

    // long logn = log2n - logDataNum;
    // long n = 1 << log2n;
    // long maskNum = logn * (logn + 1) / 2;
    // PrintUtils::printSingleMatrix("mask", mask, maskNum, n);
    // PrintUtils::printSingleMatrix("maskOther", maskOther, maskNum, n);
    // PrintUtils::printSingleMatrix("maskTable", maskTable, maskNum, n);
    // PrintUtils::printSingleMatrix("maskTableOther", maskTableOther, maskNum, n);

    plainSortingTableRecursion(ca, log2n - logDataNum, 0, 0);
}

long PlainSort::plainSortingTableRecursion(CyclicArray& ca, long logNum, long logDist, long loc) {
    if (logNum == 1) {
        compAndSwapTable(ca, loc, 1 << (logDist + logDataNum));
    } else {
        if (logDist == 0) {
            loc = plainSortingTableRecursion(ca, logNum - 1, logDist, loc);
        }
        loc = plainSortingTableRecursion(ca, logNum - 1, logDist + 1, loc);
        compAndSwapTable(ca, loc, 1 << (logDist + logDataNum));
    }
    return loc + 1;
}

void PlainSort::compAndSwapTable(CyclicArray& ca, long loc, long dist, bool _increase) {
    increase = _increase;
    setInc(increase);
    compAndSwapTable(ca, loc, dist);
}

void PlainSort::compAndSwapTable(CyclicArray& ca, long loc, long dist) {
    long length = ca.length;
    CyclicArray maskPoly(mask[inc][loc], length);
    CyclicArray maskRightPoly(maskRight[inc][loc], length);
    CyclicArray maskTablePoly(maskTable[inc][loc], length);
    CyclicArray maskTableRightPoly(maskTableRight[inc][loc], length);
    CyclicArray ca1, ca1Right, caTable, caTableRight;
    mult(ca1, ca, maskPoly);
    mult(ca1Right, ca, maskRightPoly);
    mult(caTable, ca, maskTablePoly);
    mult(caTableRight, ca, maskTableRightPoly);

    // caTable.printAsVector();

    ca.sub(ca1);
    ca.sub(ca1Right);

    if(increase) {
        ca1.rightRotate(dist);
        caTable.rightRotate(dist);
    } 
    else {
        ca1.leftRotate(dist);
        caTable.leftRotate(dist);
    }

    // cout << "before comp" << endl;
    // caTable.printAsVector();
    // caTableRight.printAsVector();
    
    minMaxTable(ca1, ca1Right, caTable, caTableRight, logDataNum, colNum, maskRight[inc][loc], maskTableRight[inc][loc]);

    if(increase) ca1.leftRotate(dist);
    else ca1.rightRotate(dist);

    ca.add(ca1);
    ca.add(ca1Right);
}

void PlainSort::minMaxTable(CyclicArray& caLeft, CyclicArray& caRight, CyclicArray& caTableLeft, CyclicArray& caTableRight, long logDataNum, long colNum, double* maskRight, double* maskTableRight) {
    long length = caLeft.length;
    for (int i = 0; i < length; i++) {
        if(caTableLeft.get(i) >= caTableRight.get(i)) {
            caTableLeft.set(i, 1);
            caTableRight.set(i, 0);
        } else {
            caTableLeft.set(i, 0);
            caTableRight.set(i, 1);
        }
    }
    CyclicArray maskTablePoly(maskTableRight, length);
    caTableLeft.mult(maskTablePoly);
    caTableRight.mult(maskTablePoly);

    caTableLeft.leftRotate(colNum);
    caTableRight.leftRotate(colNum);

    for (int i = 0; i < logDataNum; i++) {
        CyclicArray tmpMinTable(caTableLeft);
        tmpMinTable.rightRotate(1 << i);
        CyclicArray tmpMaxTable(caTableRight);
        tmpMaxTable.rightRotate(1 << i);
        caTableLeft.add(tmpMinTable);
        caTableRight.add(tmpMaxTable);
    }

    CyclicArray caTableLeftFlip(caTableLeft);
    CyclicArray caTableRightFlip(caTableRight);
    for (int i = 0; i < length; i++) {
        caTableLeftFlip.set(i, 1 - caTableLeftFlip.get(i));
        caTableRightFlip.set(i, 1 - caTableRightFlip.get(i));
    }

    CyclicArray caLeftSmall, caLeftBig, caRightSmall, caRightBig;
    mult(caLeftSmall, caLeft, caTableLeftFlip);
    mult(caLeftBig, caLeft, caTableLeft);
    mult(caRightSmall, caRight, caTableRightFlip);
    mult(caRightBig, caRight, caTableRight);

    add(caLeft, caLeftSmall, caRightSmall);
    add(caRight, caLeftBig, caRightBig);
    
}

void PlainSort::bitonicTableMerge(CyclicArray* ca, long log2n, long logNum, long logDataNum, long colNum) {
     MaskingGenerator mg(log2n, logDataNum, colNum, true);
    double** maskCol = mg.getColNumMasking();

    MaskingGenerator mg1(log2n, logDataNum, true);
    double** mask = mg1.getBitonicMergeMasking();
    MaskingGenerator mg12(log2n, logDataNum, true);
    double** maskOther = mg12.getBitonicMergeMaskingOther();
    MaskingGenerator mgTable(log2n, logDataNum, colNum, true);
    double** maskTable = mgTable.getBitonicMergeMasking();
    MaskingGenerator mgTable2(log2n, logDataNum, colNum, true);
    double** maskTableOther = mgTable2.getBitonicMergeMaskingOther();

    long maskNum = log2n - logDataNum;
    double*** maskInc = new double**[4];
    maskInc[0] = mask;
    maskInc[1] = maskOther;
    maskInc[2] = maskTable;
    maskInc[3] = maskTableOther;
    for(int i = 0; i < 4; i++) {
        mg.printMask(maskInc[i], maskNum);
    }

    MaskingGenerator mg2(log2n, logDataNum, false);
    double** maskDec1 = mg2.getBitonicMergeMasking();
    MaskingGenerator mg22(log2n, logDataNum, false);
    double** maskOtherDec = mg22.getBitonicMergeMaskingOther();
    MaskingGenerator mgTable22(log2n, logDataNum, colNum, false);
    double** maskTableDec = mgTable22.getBitonicMergeMasking();
    MaskingGenerator mgTable222(log2n, logDataNum, colNum, false);
    double** maskTableOtherDec = mgTable222.getBitonicMergeMaskingOther();
    double*** maskDec = new double**[4];
    maskDec[0] = maskDec1;
    maskDec[1] = maskOtherDec;
    maskDec[2] = maskTableDec;
    maskDec[3] = maskTableOtherDec;
    for(int i = 0; i < 4; i++) {
        mg.printMask(maskDec[i], maskNum);
    }

    bitonicTableMergeRec(ca, log2n, 0, logNum, logDataNum, colNum, maskCol, maskInc, maskDec, true);
}

void PlainSort::bitonicTableMergeRec(CyclicArray* ca, long log2n, long start, long logNum, long logDataNum, long colNum, double** maskCol, double*** maskInc, double*** maskDec, bool increase) {
    if (logNum == 0) return;        
    
    // cout << "logNum = " << logNum << endl;

    bitonicTableMergeRec(ca, log2n, start, logNum - 1, logDataNum, colNum, maskCol, maskInc, maskDec, true);
    bitonicTableMergeRec(ca, log2n, start + (1 << (logNum - 1)), logNum - 1, logDataNum, colNum, maskCol, maskInc, maskDec, false);
    
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
                CyclicArray maskminMaxTablePoly(maskCol[1], ca[0].length);
                CyclicArray caLeftTable(ca[start + left]);
                CyclicArray caRightTable(ca[start + right]);
                caLeftTable.mult(maskminMaxTablePoly);
                caRightTable.mult(maskminMaxTablePoly);
                
                minMaxTable(ca[start + left], ca[start + right], caLeftTable, caRightTable, logDataNum, colNum, maskCol[0], maskCol[1]);
            }
        }
    }
    
    for(int i = 0; i < (1 << logNum); i++) {
        // cout << "self " << start + i << ", " << increase << endl;
        double*** mask;
        if (increase) mask = maskInc;
        else          mask = maskDec;
        
        selfBitonicTableMerge(ca[start + i], log2n, logDataNum, colNum, mask, increase);
    } 
    // cout << " -- end " << logNum << " -- " << endl;
}

void PlainSort::selfBitonicTableMerge(CyclicArray& ca, long log2n, long logDataNum, long colNum, double*** mask, bool increase) {
    for(int i = 0; i <log2n - logDataNum; i++) {
        compAndSwapTable(ca, i, 1 << (log2n - 1 - i));
    }
}