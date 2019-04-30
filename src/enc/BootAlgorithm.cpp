#include "BootAlgorithm.h"

BootAlgo::BootAlgo(Parameter _param, long _sqrtIter, bool _increase) {
    param = _param;
    sqrtIter = _sqrtIter;
    invIter = sqrtIter;
    compIter = invIter;
    increase = _increase;
}

BootAlgo::BootAlgo(Parameter _param, long _invIter, long _compIter, bool _increase) {
    param = _param;
    sqrtIter = _invIter;
    invIter = _invIter;
    compIter = _compIter;
    increase = _increase;
}

BootAlgo::BootAlgo(Parameter _param, long _sqrtIter, long _invIter, long _compIter, bool _increase) {
    param = _param;
    sqrtIter = _sqrtIter;
    invIter = _invIter;
    compIter = _compIter;
    increase = _increase;
}

void BootAlgo::approxSqrt(Ciphertext& cipher, BootScheme& scheme, BootHelper& bootHelper) {
    PrintUtils::nprint("start BootAlgo::sqrt", WANT_TO_PRINT);

    Ciphertext b = cipher;
    long logp = param.logp;    
    scheme.addConstAndEqual(b, -1.0, logp);

    // one sqrtIteration : 2logp 
    // total = 2 * d * logp
    Ciphertext dummy;
    for(int i = 0; i < sqrtIter; i++) {
        PrintUtils::nprint(to_string(i) + "/" + to_string(sqrtIter - 1) + "th sqrtIteration", WANT_TO_PRINT);

        scheme.checkAndBoot(cipher, cipher.logq - param.logp < param.logq, bootHelper, param);
        scheme.checkAndBoot(b, b.logq - param.logp < param.logq, bootHelper, param);


        // make dummy = 1 - b / 2
        dummy = scheme.divByPo2(b, 1); // b - 1
        scheme.negateAndEqual(dummy); 
        scheme.addConstAndEqual(dummy, 1.0, logp); // dummy - 1

        // Update a
        // a <- a * (1 - b / 2)
        scheme.modDownToAndEqualModified(cipher, dummy, bootHelper, param); // a - 1
        scheme.multAndEqualWithBoot(cipher, dummy, bootHelper, param);
        scheme.reScaleByAndEqual(cipher, logp); // a - logp + 1

        // make dummy = (b - 3) / 4
        dummy = scheme.addConst(b, -3.0, logp);
        scheme.divByPo2AndEqual(dummy, 2); // dummy - 3

        //update b
        // b<- b * b * (b - 3) / 4
        scheme.squareAndEuqalWithBoot(b, bootHelper, param);
        scheme.reScaleByAndEqual(b, logp); // b - logp
        scheme.modDownToAndEqualModified(dummy, b, bootHelper, param);
        scheme.multAndEqualWithBoot(b, dummy, bootHelper, param);
        scheme.reScaleByAndEqual(b, logp); // b - 2logp
        scheme.modDownToAndEqualModified(cipher, b, bootHelper, param);

        PrintUtils::nprint("After sqrtIter : logq = " + to_string(cipher.logq), WANT_TO_PRINT);
    }
    PrintUtils::nprint("end Sqrt, logq = " + to_string(cipher.logq), WANT_TO_PRINT);
}

void BootAlgo::approxInverse(Ciphertext& cipher, BootScheme& scheme, BootHelper& bootHelper) {
    scheme.negateAndEqual(cipher);
    Ciphertext a = scheme.addConst(cipher, 2.0, param.logp);
    Ciphertext b = scheme.addConst(cipher, 1.0, param.logp);

    for (int _ = 0; _ < invIter; _++) {
        scheme.squareAndEuqalWithBoot(b, bootHelper, param);
        scheme.reScaleByAndEqual(b, param.logp);
        Ciphertext bPlusOne = scheme.addConst(b, 1.0, param.logp);
        scheme.modDownToAndEqualModified(a, bPlusOne, bootHelper, param);
        scheme.multAndEqualWithBoot(a, bPlusOne, bootHelper, param);
        scheme.reScaleByAndEqual(a, param.logp);
        scheme.modDownToAndEqualModified(b, a, bootHelper, param);
    }
    cipher = a;
}

void BootAlgo::minMax(Ciphertext& minCipher, Ciphertext& maxCipher, BootScheme& scheme, BootHelper& bootHelper) {
    PrintUtils::nprint("start minMax with logq = " + to_string(minCipher.logq) + ", " + to_string(maxCipher.logq), WANT_TO_PRINT);
    Ciphertext x = scheme.add(minCipher, maxCipher);
    Ciphertext y = scheme.sub(minCipher, maxCipher);
    scheme.divByPo2AndEqual(x, 1); // x - logp + 1
    scheme.divByPo2AndEqual(y, 1); // y - logp + 1
    
    scheme.squareAndEuqalWithBoot(y, bootHelper, param);
    scheme.reScaleByAndEqual(y, param.logp); // y - logp + 1

    // sqrtCipher - (2 * sqrtIter + 1) * logp + 1
    approxSqrt(y, scheme, bootHelper);

    // scheme.modDownToAndEqual(x, sqrtCipher.logq);
    scheme.modDownToAndEqualModified(x, y, bootHelper, param);

    maxCipher = scheme.add(x, y);
    minCipher = scheme.sub(x, y);
    PrintUtils::nprint("end minMax", WANT_TO_PRINT);
}

void BootAlgo::comparison(Ciphertext& cipher1, Ciphertext& cipher2, BootScheme& scheme, BootHelper& bootHelper) {
    Ciphertext sum = scheme.add(cipher1, cipher2);
    scheme.divByPo2AndEqual(sum, 1);
    approxInverse(sum, scheme, bootHelper);
    
    scheme.divByPo2AndEqual(cipher1, 1);
    scheme.modDownToAndEqualModified(cipher1, sum, bootHelper, param);
    scheme.multAndEqualWithBoot(cipher1, sum, bootHelper, param);
    scheme.reScaleByAndEqual(cipher1, param.logp);

    Ciphertext a = cipher1;
    Ciphertext b = scheme.negate(a);
    scheme.addConstAndEqual(b, 1.0, param.logp);

    for (int _ = 0; _ < compIter; _++) {
        for (int i = 0; i < 2; i++) {
            scheme.squareAndEuqalWithBoot(a, bootHelper, param);
            scheme.reScaleByAndEqual(a, param.logp);
            scheme.squareAndEuqalWithBoot(b, bootHelper, param);
            scheme.reScaleByAndEqual(b, param.logp);
        }
        
        Ciphertext inv = scheme.add(a, b);
        approxInverse(inv, scheme, bootHelper);
        scheme.modDownToAndEqualModified(a, inv, bootHelper, param);
        scheme.multAndEqualWithBoot(a, inv, bootHelper, param);
        scheme.reScaleByAndEqual(a, param.logp);

        // scheme.modDownToAndEqualModified(b, a, bootHelper, param);
        b = scheme.negate(a);
        scheme.addConstAndEqual(b, 1.0, param.logp);
    }
    cipher1 = a;
    cipher2 = b;
}

void BootAlgo::compAndSwap(Ciphertext& cipher, double* mask, long dist, BootScheme& scheme, Ring& ring, BootHelper& bootHelper) {
    PrintUtils::nprint("start compAndSwap with logq = " + to_string(cipher.logq), WANT_TO_PRINT);
    ZZ* maskPoly = new ZZ[1 << param.logN];
    ring.encode(maskPoly, mask, cipher.n, param.logp);
    Ciphertext dummy = cipher;
    scheme.multByPolyAndEqualWithBoot(dummy, maskPoly, bootHelper, param);
    scheme.reScaleByAndEqual(dummy, param.logp);
    scheme.modDownToAndEqualModified(cipher, dummy, bootHelper, param);
    scheme.subAndEqual(cipher, dummy);
    if(increase) scheme.rightRotateFastAndEqual(dummy, dist);
    else         scheme.leftRotateFastAndEqual(dummy, dist); 
    
    minMax(dummy, cipher, scheme,  bootHelper);
    
    if(increase) scheme.leftRotateFastAndEqual(dummy, dist); 
    else         scheme.rightRotateFastAndEqual(dummy, dist);
        
    scheme.addAndEqual(cipher, dummy);   

    PrintUtils::nprint("end compAndSwap", WANT_TO_PRINT);
}

void BootAlgo::selfBitonicMerge(Ciphertext& cipher, double** mask, BootScheme& scheme, Ring& ring, BootHelper& bootHelper) {
    for(int i = 0; i < param.log2n; i++) {
        compAndSwap(cipher, mask[i], 1 << (param.log2n - 1 - i), scheme, ring, bootHelper);
    }
}

void BootAlgo::reverse(Ciphertext& cipher, double** mask, BootScheme& scheme, Ring& ring, BootHelper& bootHelper) {
    scheme.leftRotateFastAndEqual(cipher, 1 << (param.log2n - 1));
    for(int i = 1; i < param.log2n; i++) {
        ZZ* maskPoly = new ZZ[1 << param.logN];
        ring.encode(maskPoly, mask[i], cipher.n, param.logp);
        Ciphertext dummy = cipher;
        scheme.multByPolyAndEqualWithBoot(dummy, maskPoly, bootHelper, param);
        scheme.reScaleByAndEqual(dummy, param.logp);
        scheme.modDownToAndEqualModified(cipher, dummy, bootHelper, param);
        scheme.subAndEqual(cipher, dummy);

        scheme.leftRotateFastAndEqual(cipher, 1 << (param.log2n - i - 1));
        scheme.rightRotateFastAndEqual(dummy, 1 << (param.log2n - i - 1));
        scheme.addAndEqual(cipher, dummy);
    }
}

void BootAlgo::compAndSwapTable(Ciphertext& cipher, long logDataNum, double* mask, double* maskRight, double* maskTable, double* maskTableRight, long dist, BootScheme& scheme, Ring& ring, BootHelper& bootHelper) {
    std::cout <<  "run compAndSwapTable" << endl;
    ZZ* maskPoly = new ZZ[1 << param.logN]();
    ZZ* maskRightPoly = new ZZ[1 << param.logN]();
    ZZ* maskTablePoly = new ZZ[1 << param.logN]();
    ZZ* maskTableRightPoly = new ZZ[1 << param.logN]();
    std::cout << "check1" << std::endl;
    ring.encode(maskPoly,           mask,           cipher.n, param.logp);
    std::cout << "check2" << std::endl;
    ring.encode(maskRightPoly,      maskRight,      cipher.n, param.logp);
    std::cout << "check3" << std::endl;
    ring.encode(maskTablePoly,      maskTable,      cipher.n, param.logp);
    std::cout << "check4" << std::endl;
    ring.encode(maskTableRightPoly, maskTableRight, cipher.n, param.logp);
    std::cout << "check5" << std::endl;
    Ciphertext cipherTable = cipher;
    Ciphertext cipherTableRight = cipher;
    scheme.multByPolyAndEqualWithBoot(cipherTable, maskTablePoly, bootHelper, param);
    scheme.multByPolyAndEqualWithBoot(cipherTableRight, maskTableRightPoly, bootHelper, param);
    scheme.reScaleByAndEqual(cipherTable, param.logp);
    scheme.reScaleByAndEqual(cipherTableRight, param.logp);

    // scheme.decryptAndPrint("table", secretKey, cipherTable);
    // scheme.decryptAndPrint("tableRight", secretKey, cipherTableRight);
    std::cout << "check" << std::endl;
    Ciphertext cipher1 = cipher;
    Ciphertext cipher1Right = cipher;
    scheme.multByPolyAndEqualWithBoot(cipher1, maskPoly, bootHelper, param);
    scheme.multByPolyAndEqualWithBoot(cipher1Right, maskRightPoly, bootHelper, param);
    scheme.reScaleByAndEqual(cipher1, param.logp);
    scheme.reScaleByAndEqual(cipher1Right, param.logp);
    std::cout << "check" << std::endl;
    // scheme.decryptAndPrint("1", secretKey, cipher1);
    // scheme.decryptAndPrint("1Right", secretKey, cipher1Right);
    
    scheme.modDownToAndEqualModified(cipher, cipher1, bootHelper, param);
    scheme.subAndEqual(cipher, cipher1);
    scheme.subAndEqual(cipher, cipher1Right);

    // scheme.decryptAndPrint("cipher", secretKey, cipher);

    scheme.rightRotateFastAndEqual(cipherTable, dist);
    
    comparison(cipherTable, cipherTableRight, scheme, bootHelper);
    
    scheme.multByPolyAndEqualWithBoot(cipherTableRight, maskTableRightPoly, bootHelper, param);
    scheme.reScaleByAndEqual(cipherTableRight, param.logp);
    scheme.modDownToAndEqualModified(cipherTable, cipherTableRight, bootHelper, param);
    
    // scheme.decryptAndPrint("table", secretKey, cipherTable);
    // scheme.decryptAndPrint("tableRight", secretKey, cipherTableRight);
    
    scheme.leftRotateFastAndEqual(cipherTable, dist);    

    

    // TODO: If colNum changed -> this will occur problems
    for (int i = 0; i < logDataNum; i++) {
        Ciphertext tmpTable = scheme.rightRotateFast(cipherTable, 1 << i);
        Ciphertext tmpTableRight = scheme.rightRotateFast(cipherTableRight, 1 << i);
        scheme.addAndEqual(cipherTable, tmpTable);
        scheme.addAndEqual(cipherTableRight, tmpTableRight);
    }

    // cout << "after rot" << endl;
    // scheme.decryptAndPrint("table", secretKey, cipherTable);
    // scheme.decryptAndPrint("tableRight", secretKey, cipherTableRight);

    Ciphertext cipherTableFlip = scheme.negate(cipherTable);
    Ciphertext cipherTableFlipRight = scheme.negate(cipherTableRight);

    scheme.addByPolyAndEqual(cipherTableFlip, maskPoly, param.logp);
    scheme.addByPolyAndEqual(cipherTableFlipRight, maskRightPoly, param.logp);

    // scheme.decryptAndPrint("tableFlip", secretKey, cipherTableFlip);
    // scheme.decryptAndPrint("tableFlipRight", secretKey, cipherTableFlipRight);

    scheme.modDownToAndEqualModified(cipher1, cipherTable, bootHelper, param);
    scheme.modDownToAndEqualModified(cipher1Right, cipherTableRight, bootHelper, param);
    
    Ciphertext cipherLeftSmall = cipher1;
    scheme.multAndEqualWithBoot(cipherLeftSmall, cipherTableFlip, bootHelper, param);
    scheme.reScaleByAndEqual(cipherLeftSmall, param.logp);

    Ciphertext cipherLeftBig = cipher1;
    scheme.multAndEqualWithBoot(cipherLeftBig, cipherTable, bootHelper, param);
    scheme.reScaleByAndEqual(cipherLeftBig, param.logp);

    Ciphertext cipherRightSmall = cipher1Right;
    scheme.multAndEqualWithBoot(cipherRightSmall, cipherTableFlipRight, bootHelper, param);
    scheme.reScaleByAndEqual(cipherRightSmall, param.logp);

    Ciphertext cipherRightBig = cipher1Right;
    scheme.multAndEqualWithBoot(cipherRightBig, cipherTableRight, bootHelper, param);
    scheme.reScaleByAndEqual(cipherRightBig, param.logp);
 
    // scheme.decryptAndPrint("LeftSmall", secretKey, cipherLeftSmall);
    // scheme.decryptAndPrint("LeftBig", secretKey, cipherLeftBig);
    // scheme.decryptAndPrint("RightSmall", secretKey, cipherRightSmall);
    // scheme.decryptAndPrint("RihtBig", secretKey, cipherRightBig);

    scheme.rightRotateAndEqual(cipherLeftBig, dist);
    scheme.leftRotateAndEqual(cipherRightSmall, dist);
    scheme.modDownToAndEqualModified(cipher, cipherLeftSmall, bootHelper, param);

    scheme.addAndEqual(cipher, cipherLeftSmall);
    scheme.addAndEqual(cipher, cipherLeftBig);
    scheme.addAndEqual(cipher, cipherRightSmall);
    scheme.addAndEqual(cipher, cipherRightBig);
}