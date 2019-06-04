#include "EncSorting.h"

EncSorting::EncSorting(Parameter _param, long* iter, long numOfIter, bool _increase, bool _printCondition) : param(_param), increase(_increase) {
    encAlgo = new EncAlgo(_param, iter, numOfIter, _printCondition);
    setInc(_increase);
    if(numOfIter > 2) {
        logDataNum = iter[0];
        colNum = iter[1];
        plainSort = new PlainSort(param.log2n, logDataNum, colNum);
    } else {
        plainSort = new PlainSort(param.log2n);
    }
}

void EncSorting::showDiffFromPlain() {
    showDiff = true;
}

Ciphertext EncSorting::encrypt(double* mvec) {
    ca = new CyclicArray(mvec, 1 << param.log2n);
    return encAlgo->encrypt(mvec);
}
complex<double>* EncSorting::decrypt(Ciphertext& cipher) {
    return encAlgo->decrypt(cipher);
}

void EncSorting::genMaskPoly() {
    if(maskPolyGen) return;
    maskPolyGen = true;

    maskPoly = new ZZ**[2];
    for (int inc = 0; inc < 2; inc++) {
        MaskingGenerator mg(param.log2n, inc == 0);
        double** mask = mg.getMasking();
        maskPoly[inc] = new ZZ*[mg.maskNum];
        for (int i = 0; i < mg.maskNum; i++) {
            maskPoly[inc][i] = encAlgo->encode(mask[i]);
        }    
    }
}

void EncSorting::genTableMaskPoly() {
    if(maskTableGen) return;
    maskTableGen = true;

    maskPoly = new ZZ**[2];
    for (int inc = 0; inc < 2; inc++) {
        MaskingGenerator mg(param.log2n, logDataNum, inc == 0);
        double** mask = mg.getMasking();
        maskPoly[inc] = new ZZ*[mg.maskNum];
        for (int i = 0; i < mg.maskNum; i++) {
            maskPoly[inc][i] = encAlgo->encode(mask[i]);
        }    
    }
    
    maskRightPoly = new ZZ**[2];
    for (int inc = 0; inc < 2; inc++) {
        MaskingGenerator mg(param.log2n, logDataNum, inc == 0);
        double** mask = mg.getMaskingOther();
        maskRightPoly[inc] = new ZZ*[mg.maskNum];
        for (int i = 0; i < mg.maskNum; i++) {
            maskRightPoly[inc][i] = encAlgo->encode(mask[i]);
        }    
    }

    maskTablePoly = new ZZ**[2];
    for (int inc = 0; inc < 2; inc++) {
        MaskingGenerator mg(param.log2n, logDataNum, colNum, inc == 0);
        double** mask = mg.getMasking();
        maskTablePoly[inc] = new ZZ*[mg.maskNum];
        for (int i = 0; i < mg.maskNum; i++) {
            maskTablePoly[inc][i] = encAlgo->encode(mask[i]);
        }    
    }

    maskTableRightPoly = new ZZ**[2];
    for (int inc = 0; inc < 2; inc++) {
        MaskingGenerator mg(param.log2n, logDataNum, colNum, inc == 0);
        double** mask = mg.getMaskingOther();
        maskTableRightPoly[inc] = new ZZ*[mg.maskNum];
        for (int i = 0; i < mg.maskNum; i++) {
            maskTableRightPoly[inc][i] = encAlgo->encode(mask[i]);
        }    
    }
}

void EncSorting::genMergeMaskingPoly() {
    if(maskMergeGen) return;
    maskMergeGen = true;

    maskMergePoly = new ZZ**[2];
    for (int inc = 0; inc < 2; inc++) {
        MaskingGenerator mg(param.log2n, inc == 0);
        double** mask = mg.getBitonicMergeMasking();
        maskMergePoly[inc] = new ZZ*[mg.maskNum];
        for (int i = 0; i < mg.maskNum; i++) {
            maskMergePoly[inc][i] = encAlgo->encode(mask[i]);
        }    
    }
}

void EncSorting::genTableMergeMaskingPoly() {
    if(maskMergeGen) return;
    maskMergeGen = true;

    maskMergePoly = new ZZ**[2];
    for (int inc = 0; inc < 2; inc++) {
        MaskingGenerator mg(param.log2n, logDataNum, inc == 0);
        double** mask = mg.getBitonicMergeMasking();
        maskMergePoly[inc] = new ZZ*[mg.maskNum];
        for (int i = 0; i < mg.maskNum; i++) {
            maskMergePoly[inc][i] = encAlgo->encode(mask[i]);
        }    
    }

    maskMergeColPoly = new ZZ**[2];
    for (int inc = 0; inc < 2; inc++) {
        MaskingGenerator mg(param.log2n, logDataNum, colNum, inc == 0);
        double** mask = mg.getBitonicMergeMasking();
        maskMergeColPoly[inc] = new ZZ*[mg.maskNum];
        for (int i = 0; i < mg.maskNum; i++) {
            maskMergeColPoly[inc][i] = encAlgo->encode(mask[i]);
        }    
    }

    maskMergeRightPoly = new ZZ**[2];
    for (int inc = 0; inc < 2; inc++) {
        MaskingGenerator mg(param.log2n, logDataNum, inc == 0);
        double** mask = mg.getBitonicMergeMasking();
        maskMergeRightPoly[inc] = new ZZ*[mg.maskNum];
        for (int i = 0; i < mg.maskNum; i++) {
            maskMergeRightPoly[inc][i] = encAlgo->encode(mask[i]);
        }    
    }

    maskMergeRightColPoly = new ZZ**[2];
    for (int inc = 0; inc < 2; inc++) {
        MaskingGenerator mg(param.log2n, logDataNum, inc == 0);
        double** mask = mg.getBitonicMergeMasking();
        maskMergeRightColPoly[inc] = new ZZ*[mg.maskNum];
        for (int i = 0; i < mg.maskNum; i++) {
            maskMergeRightColPoly[inc][i] = encAlgo->encode(mask[i]);
        }    
    }
}

void EncSorting::genMinMaxTablePoly() {
    if(maskMinMaxTableGen) return;
    maskMinMaxTableGen = true;

    maskMinMaxTablePoly = new ZZ*[2];
    MaskingGenerator mg(param.log2n, logDataNum, colNum, true);
    double** mask = mg.getColNumMasking();
    for (int i = 0; i < 2; i++) {
        maskMinMaxTablePoly[i] = new ZZ[1 << param.logN];
        maskMinMaxTablePoly[i] = encAlgo->encode(mask[i]);
    }
}

void EncSorting::setInc(bool _increase) {
    increase = _increase;
    if (increase) {
        inc = 0;
    } else {
        inc = 1;
    }
}

void EncSorting::runEncSorting(Ciphertext& cipher, bool _increase) {
    increase = _increase;
    setInc(increase);
    runEncSorting(cipher);
}

void EncSorting::runEncSorting(Ciphertext& cipher) {
    genMaskPoly();
    sortingRecursion(cipher, param.log2n, 0, 0);
    encAlgo->showTotalCount();
}

long EncSorting::sortingRecursion(Ciphertext& cipher, long logNum, long logDist, long loc) {
    TimeUtils timeutils;
    
    PrintUtils::nprint("run encBatcherSort with loc = " + to_string(loc), WANT_TO_PRINT);
    if (logNum == 1) {
        timeutils.start(to_string(loc)+"th CompAndSwap");
        encSwapAndCompare(cipher, logDist, loc);
        timeutils.stop(to_string(loc)+"th CompAndSwap");
    } else {
        if (logDist == 0) {
            loc = sortingRecursion(cipher, logNum - 1, logDist, loc);
        }
        loc = sortingRecursion(cipher, logNum - 1, logDist + 1, loc);
        timeutils.start(to_string(loc)+"th CompAndSwap");
        encSwapAndCompare(cipher, logDist, loc);
        timeutils.stop(to_string(loc)+"th CompAndSwap");
    }
    return loc + 1;
}

void EncSorting::encSwapAndCompare(Ciphertext& cipher, long logDist, long loc) {
    cout << "start encSwapAndCompare" << endl;
    cout << "start encSwap" << endl;
    encAlgo->encSwap(cipher, maskPoly[inc][loc], 1 << logDist, increase);
    cout << "end encSwap" << endl;

    if (showDiff) {
        cout << "start compAndSwap" << endl;
        plainSort->compAndSwap(*ca, loc, 1 << logDist, increase);
        cout << "end compAndSwap" << endl;
        double* plain = ca->getArray();
        complex<double>* enc = encAlgo->decrypt(cipher);
        cout << loc << "th compAndSwap Result" << endl;
        PrintUtils::printFewArrays(plain, enc, cipher.n);
        cout << endl << "******************" << endl;
        PrintUtils::averageDifference(plain, enc, cipher.n);
        cout << "******************" << endl << endl;
    }
    encAlgo->showCurrentCountAndReset();
}

void EncSorting::runEncTableSorting(Ciphertext& cipher) {
    genTableMaskPoly();
    sortingTableRecursion(cipher, param.log2n - logDataNum, 0, 0);
    encAlgo->showTotalCount();
}

long EncSorting::sortingTableRecursion(Ciphertext& cipher, long logNum, long logDist, long loc) {
    TimeUtils timeutils;
    if (logNum == 1) {
        timeutils.start(to_string(loc)+"th CompAndSwap");
        encSwapTableAndCompare(cipher, 1 << (logDist + logDataNum), loc); 
        timeutils.stop(to_string(loc)+"th CompAndSwap with " + to_string(logNum) + ", " + to_string(logDist));
    } else {
        if (logDist == 0) {
            loc = sortingTableRecursion(cipher, logNum - 1, logDist, loc);
        }
        loc = sortingTableRecursion(cipher, logNum - 1, logDist + 1, loc);
        timeutils.start(to_string(loc)+"th CompAndSwap");
        encSwapTableAndCompare(cipher, 1 << (logDist + logDataNum), loc); 
        timeutils.stop(to_string(loc)+"th CompAndSwap with " + to_string(logNum) + ", " + to_string(logDist));
    }
    return loc + 1;
}

void EncSorting::encSwapTableAndCompare(Ciphertext& cipher, long dist, long loc) {
    encAlgo->encSwapTable(cipher, maskPoly[inc][loc], maskRightPoly[inc][loc], maskTablePoly[inc][loc], maskTableRightPoly[inc][loc], dist, increase);

    if (showDiff) {
        cout << "start compAndSwap" << endl;
        plainSort->compAndSwapTable(*ca, loc, dist, increase);
        cout << "end compAndSwap" << endl;
        double* plain = ca->getArray();
        complex<double>* enc = encAlgo->decrypt(cipher);
        cout << loc << "th compAndSwap Result" << endl;
        PrintUtils::printFewArrays(plain, enc, cipher.n);
        cout << endl << "******************" << endl;
        PrintUtils::averageDifference(plain, enc, cipher.n);
        cout << "******************" << endl << endl;
    }
    encAlgo->showCurrentCountAndReset();
}

void EncSorting::bitonicMerge(Ciphertext* cipher, long logNum) {
    genMergeMaskingPoly();
    bitonicMergeRec(cipher, 0, logNum, true);
    encAlgo->showTotalCount();
}

void EncSorting::bitonicMergeRec(Ciphertext* cipher, long start, long logNum, bool increase) {
    if (logNum == 0) return;        

    bitonicMergeRec(cipher, start, logNum - 1, true);
    bitonicMergeRec(cipher, start + (1 << (logNum - 1)), logNum - 1, false);

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
                cout << "-- Run minMax (" << start + left << ", " << start + right << ")" << endl;
                encAlgo->minMaxAlgorithm(cipher[start + left], cipher[start + right]);
            }
        }
    }
    encAlgo->showCurrentCountAndReset();
        cout << "----- run selfBitonicMerge " << endl;
    for(int i = 0; i < (1 << logNum); i++) {        
        cout << "           -- ctxt " << start + i << endl;
        // encAlgo->selfBitonicMerge(cipher[start + i], increase);
    } 
    encAlgo->showCurrentCountAndReset();
    cout << "-----";
}


void EncSorting::bitonicTableMerge(Ciphertext* cipher, long logNum) {
    genMinMaxTablePoly();
    genTableMergeMaskingPoly();
    // MaskingGenerator mg(param.log2n, logDataNum, colNum, true);
    // double** maskCol = mg.getColNumMasking();
    // ZZ** maskminMaxTablePoly = new ZZ*[2];
    // maskminMaxTablePoly[0] = new ZZ[1 << param.logN];
    // maskminMaxTablePoly[1] = new ZZ[1 << param.logN];
    // ring.encode(maskminMaxTablePoly[0], maskCol[0], 1 << param.log2n, param.logp);
    // ring.encode(maskminMaxTablePoly[1], maskCol[1], 1 << param.log2n, param.logp);


    // MaskingGenerator mg1(param.log2n, logDataNum, true);
    // double** mask = mg1.getBitonicMergeMasking();
    // MaskingGenerator mg12(param.log2n, logDataNum, true);
    // double** maskOther = mg12.getBitonicMergeMaskingOther();
    // MaskingGenerator mgTable(param.log2n, logDataNum, colNum, true) {
    // double** maskTableOther = mgTable2.getBitonicMergeMaskingOther();

    // long maskNum = param.log2n - logDataNum;
    // double*** maskInc = new double**[4];
    // maskInc[0] = mask;
    // maskInc[1] = maskOther {
    // // for(int i = 0; i < 4; i++) {
    // //     mg.printMask(maskInc[i], maskNum);
    // // }

    // MaskingGenerator mg2(param.log2n, logDataNum, false);
    // double** maskDec1 = mg2.getBitonicMergeMasking();
    // MaskingGenerator mg22(param.log2n, logDataNum, false);
    // double** maskOtherDec = mg22.getBitonicMergeMaskingOther();
    // MaskingGenerator mgTable22(param.log2n, logDataNum, colNum, false);
    // double** maskTableDec = mgTable22.getBitonicMergeMasking();
    // MaskingGenerator mgTable222(param.log2n, logDataNum, colNum, false);
    // double** maskTableOtherDec = mgTable222.getBitonicMergeMaskingOther();
    // double*** maskDec = new double**[4];
    // maskDec[0] = maskDec1;
    // maskDec[1] = maskOtherDec;
    // maskDec[2] = maskTableDec;
    // maskDec[3] = maskTableOtherDec;
    // // for(int i = 0; i < 4; i++) {
    // //     mg.printMask(maskDec[i], maskNum);
    // // }

    // bitonicTableMergeRec(cipher, 0, logNum, logDataNum, colNum, maskminMaxTablePoly, maskInc, maskDec, scheme, ring, bootHelper, sk, true);
    bitonicTableMergeRec(cipher, 0, logNum, true);

    encAlgo-> showTotalCount();
}

void EncSorting::bitonicTableMergeRec(Ciphertext* cipher, long start, long logNum, bool increase) {
    setInc(increase);

    cout << "TableMerge " << start << ", " << logNum << endl;
    if (logNum == 0) return;        
    bitonicTableMergeRec(cipher, start, logNum - 1, true);
    bitonicTableMergeRec(cipher, start + (1 << (logNum - 1)), logNum - 1, false);

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
                // Ciphertext cipherLeftTable = encAlgo->scheme->multByPolyWithBoot(cipher[start + left], maskMinMaxTablePoly[1], *encAlgo->bootHelper, param);
                // Ciphertext cipherLeftTable = encAlgo->scheme->multByPolyWithBoot(cipher[start + right], maskMinMaxTablePoly[1], *encAlgo->bootHelper, param);
                // encAlgo->scheme->reScaleByAndEqual(cipherLeftTable, param.logp);
                // encAlgo->scheme->reScaleByAndEqual(cipherRightTable, param.logp);
                // encAlgo->minMaxTable(cipher[start + left], cipher[start + right], cipherLeftTable, cipherRightTable, maskminMaxTablePoly[0], maskminMaxTablePoly[1]);
            }
        }
    }
    encAlgo->showCurrentCountAndReset();
    
    // double*** mask;   
    // if (increase) {
    //     mask = maskInc;
    // }
    // else {
    //     mask = maskDec;
    // }

    cout << "----- run selfTableMerge " << endl;
    for(int i = 0; i < (1 << logNum); i++) {        
        cout << "           -- ctxt " << start + i << endl;
        // scheme.decryptAndPrint("startCipher", sk, cipher[start + i]);
        // encAlgo->selfTableMerge(cipher[start + i], increase);
        // scheme.decryptAndPrint("EndCipher", sk, cipher[start + i]);

        cout << "           -- end selfTableMerge" << endl;
    } 
    encAlgo->showCurrentCountAndReset();
    cout << "-----";
}

void EncSorting::reverseHalf(Ciphertext* cipher, long logNum) {
    long num = 1 << logNum;
    genMergeMaskingPoly();
    // cout << "---- Reverse odd ctxt " << endl;
    for (int i = 0; i < num / 2; i++) {
        encAlgo->reverse(cipher[i], maskMergePoly[0]);
    }
    cout << "---- modDown even ctxt " << endl;
    for (int i = 0; i < num / 2; i++) {
        encAlgo->scheme->modDownToAndEqualModified(cipher[2 * i], cipher[2 * i + 1], *encAlgo->bootHelper, encAlgo->param);
    }    
}

// void EncSorting::runBitonicSortDec(Ciphertext& cipher, BootScheme& scheme, Ring& ring, BootHelper& bootHelper, bool increase, SecretKey sk) {
//     double*** maskReverse = new double**[param.log2n];
//     MaskingGenerator mg(param.log2n, increase);
//     for (int i = 1; i < param.log2n; i++) {
//         maskReverse[i-1] = mg.getReverseMasking(i);
//         cout << "i = " << i << endl;
//         mg.printMask(maskReverse[i-1], i);
//     }

//     double*** maskReverseRight = new double**[param.log2n];
//     MaskingGenerator mg2(param.log2n, increase);
//     for (int i = 1; i < param.log2n; i++) {
//         maskReverseRight[i-1] = mg2.getReverseMaskingRight(i);
//         cout << "i = " << i << endl;
//         mg2.printMask(maskReverseRight[i-1], i);
//     }

//     for (int i = 1; i < param.log2n; i++) {
//         cout << "i = " << i << endl;
//         for(int j = 0; j < i; j++) {
//             cout << "maskReverse[" << i  - 1 << "][" << j << "] = ";  
//             for(int k = 0; k < (1 << param.log2n); k++) {
//                 cout << maskReverse[i-1][j][k] << ",";
//             }cout << endl;
//         }
//     }

//     MaskingGenerator mg3(param.log2n, increase);
//     double** mask = mg3.getBitonicMergeMasking();

//     bootAlgo = BootAlgo(param, sqrtIter, increase);
//     TimeUtils timeutils;
//     for (int i = 0; i < param.log2n; i++) {
        
//         for(int j = i; j >= 0; j--) {
//             timeutils.start(to_string(j)+", " + to_string(i) + "th halfCleaner");
//             scheme.resetImagErrorAndEqual(cipher);
//             bootAlgo.halfCleaner(cipher, mask[param.log2n - 1 - j], 1 << j, scheme, ring, bootHelper, sk);
//             timeutils.stop(to_string(j)+", " + to_string(i) + "th halfCleaner");
//             scheme.showCurrentCount();
//             scheme.resetCount();
//         }
        
//         if (i < param.log2n - 1) {
//         timeutils.start(to_string(i) + "th Reverse");
//         bootAlgo.reverse(cipher, maskReverse[i], maskReverseRight[i], i+1, scheme, ring, bootHelper);
//         timeutils.stop(to_string(i) + "th Reverse");
//         scheme.showCurrentCount();
//         scheme.resetCount();
//         }
//     }
//     scheme.showTotalCount();
// }

// // void EncSorting::runBitonicSortRecursion(Ciphertext& cipher, long logNum, double** mask, BootScheme& scheme, Ring& ring, BootHelper& bootHelper) {
    
// // }
