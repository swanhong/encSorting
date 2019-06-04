#include "EncAlgo.h"

EncAlgo::EncAlgo(Parameter _param, long* _iter, long numOfIters, bool _printCondition) :
    EncInterface(_param, _printCondition), iter(_iter) {
    long start = 0;
    std::cout << "***************************" << std::endl;
    std::cout << "Other Parameters" << std::endl;
    if(numOfIters > 2) {
        logDataNum = iter[0];
        colNum = iter[1];
        std::cout << "logDataNum = " << iter[0] << ", colNum = " << iter[1] << std::endl;
        start = 2;
    }
    std::cout << "iter = ";
    for (int i = 0; i < numOfIters - start; i++) {
        std::cout << iter[start + i];
        if(i < numOfIters - start - 1) {
            std::cout << ", ";
        }
    }
    std::cout << endl;
    std::cout << "want to print? : " << printCondition << std::endl;
    std::cout << "***************************" << std::endl;
    std::cout << std::endl;
    minMaxLoc += start;
    invLoc += start;
    compLoc += start;
}

void EncAlgo::sqrtAlgorithm(Ciphertext& cipher) {
    approxSqrt(cipher);
    // approxSqrt2(cipher);
}

void EncAlgo::minMaxAlgorithm(Ciphertext& minCipher, Ciphertext& maxCipher) {
    // minMax(minCipher, maxCipher);
    newMinMax(minCipher, maxCipher);
}

void EncAlgo::comparisonAlgorithm(Ciphertext& a, Ciphertext& b) {
    // comparison(a, b);
    newComparison(a, b);
}

void EncAlgo::evalFcn(Ciphertext& cipher) {
    
    // x <- x (3 - x^2) / 2, consumes 2 * logq
    for(int i = 0; i < iter[minMaxLoc]; i++) {
        scheme->checkModulusAndBoot(cipher, 2 * param.logp + 1, *bootHelper, param);
        scheme->resetImagErrorAndEqual(cipher);
        nprint("start cipher" + to_string(i), cipher);
        Ciphertext square = scheme->square(cipher);
        scheme->reScaleByAndEqual(square, param.logp);
        nprint("square" + to_string(i), square);
        scheme->negateAndEqual(square);
        scheme->addConstAndEqual(square, 3.0, param.logp);
        nprint("3 - x^2" + to_string(i), square);
        scheme->divByPo2AndEqual(cipher, 1);
        nprint("x / 2" + to_string(i), cipher);
        scheme->modDownToAndEqual(cipher, square.logq);
        scheme->multAndEqual(cipher, square);
        scheme->reScaleByAndEqual(cipher, param.logp);
        nprint("final x (3 - x^2) / 2" + to_string(i), cipher);
    }
    // scheme->addConstAndEqual(cipher, 1.0, param.logp);
    // scheme->divByPo2AndEqual(cipher, 1);

    // ********************

    // // x <- ( x^4 * (-5x^3 + 21x) + 35(-x^3 + x) ) / 16
    // for (int i = 0; i < sqrtIter; i++) {
    //     Ciphertext cipher2 = scheme->square(cipher);
    //     scheme->reScaleByAndEqual(cipher2, param.logp);
    //     scheme->modDownByAndEqual(cipher, param.logp);
    //     Ciphertext cipher3 = scheme->mult(cipher, cipher2);
    //     scheme->reScaleByAndEqual(cipher3, param.logp);
    //     Ciphertext cipher4 = scheme->square(cipher2);
    //     scheme->reScaleByAndEqual(cipher4, param.logp);
    //     scheme->modDownByAndEqual(cipher, param.logp);
    //     Ciphertext cipherLow = scheme->sub(cipher, cipher3);
    //     // cipherLow *= 35
    //     Ciphertext cipherLow35 = cipherLow;
    //     Ciphertext cipher_21 = cipher;
    //     Ciphertext cipher3_5 = cipher3;
    //     for (int j = 0; j < 34; j++) {
    //         scheme->addAndEqual(cipherLow35, cipherLow);
    //     }
    //     for (int j = 0; j < 20; j++) {
    //         scheme->addAndEqual(cipher_21, cipher);
    //     }
    //     for (int j = 0; j < 4; j++) {
    //         scheme->addAndEqual(cipher3_5, cipher3);
    //     }
    //     Ciphertext cipherSub = scheme->sub(cipher_21, cipher3_5);
    //     scheme->multAndEqual(cipher4, cipherSub);
    //     scheme->reScaleByAndEqual(cipher4, param.logp);
    //     scheme->modDownByAndEqual(cipherLow35, param.logp);
    //     cipher = scheme->add(cipher4, cipherLow35);
    //     scheme->divByPo2AndEqual(cipher, 4);
    // }
    // scheme->addConstAndEqual(cipher, 1.0, param.logp);
    // scheme->divByPo2AndEqual(cipher, 1);
}

void EncAlgo::approxSqrt(Ciphertext& cipher) {

    // scheme->checkAndBoot(cipher, true, *bootHelper, param);
    scheme->resetImagErrorAndEqual(cipher);

    Ciphertext b = cipher;
    long logp = param.logp;    
    scheme->addConstAndEqual(b, -1.0, logp);

    nprint("start b", b);

    // one sqrtIteration : 2logp 
    // total = 2 * d * logp
    Ciphertext dummy;
    for(int i = 0; i < iter[minMaxLoc]; i++) {
        scheme->checkModulusAndBoot(cipher, 2 * param.logp + 2, *bootHelper, param);
        scheme->checkModulusAndBoot(b, 2 * param.logp + 2, *bootHelper, param);


        // make dummy = 1 - b / 2
        dummy = scheme->divByPo2(b, 1); // b - 1
        scheme->negateAndEqual(dummy); 
        scheme->addConstAndEqual(dummy, 1.0, logp); // dummy - 1
	
        nprint(to_string(i) + "dummy = 1 - b / 2", dummy);
        
        // Update a
        // a <- a * (1 - b / 2)
        scheme->modDownToAndEqualModified(cipher, dummy, *bootHelper, param); // a - 1
        scheme->multAndEqualWithBoot(cipher, dummy, *bootHelper, param);
        scheme->reScaleByAndEqual(cipher, logp); // a - logp + 1
	
        // make dummy = (b - 3) / 4
        dummy = scheme->addConst(b, -3.0, logp);
        scheme->divByPo2AndEqual(dummy, 2); // dummy - 3
	
        nprint(to_string(i) + "dummy = (b - 3) / 4", dummy);
        
        //update b
        // b<- b * b * (b - 3) / 4
        // scheme->squareByConjugateAndEqual(b, param.logp);
        scheme->squareAndEqual(b);
        scheme->reScaleByAndEqual(b, param.logp);
        scheme->modDownToAndEqualModified(dummy, b, *bootHelper, param);
        scheme->multAndEqualWithBoot(b, dummy, *bootHelper, param);
        scheme->reScaleByAndEqual(b, logp); // b - 2logp
	
        nprint(to_string(i) + "b = b * b* (b-3)/4", b);

        if(i < iter[minMaxLoc] - 1){
            scheme->resetImagErrorAndEqual(b);
            scheme->resetImagErrorAndEqual(cipher);
            scheme->modDownToAndEqualModified(cipher, b, *bootHelper, param);
        }     

        nprint("cipher" + to_string(i), cipher);
        nprint("bModDown" + to_string(i), b);
    }
}

void EncAlgo::approxSqrt2(Ciphertext& cipher) {
    long logp = param.logp;  

    // We consumes modulus 3 * logp + 1 for each iteration. So check remaining modulus before copy cipher to r.
    scheme->checkAndBoot(cipher, cipher.logq - (2 * logp + 1) < param.logq, *bootHelper, param);
    
    Ciphertext r = cipher;
    scheme->negateAndEqual(r); 
    scheme->addAndEqual(r, cipher);
    scheme->addConstAndEqual(r, 0.5, logp);

    nprint("Start r", r);
       
    for(int i = 0; i < iter[minMaxLoc]; i++) {
	    scheme->checkAndBoot(r, r.logq - 2 * logp - 1 < param.logq, *bootHelper, param);
        scheme->checkAndBoot(cipher, cipher.logq < r.logq, *bootHelper, param); // once cipher has been bootstrpped, it never bootstrapp again while the iterations end.
        
        // w <- cipher * r^3 = (cipher * r) * r^2
        Ciphertext w = scheme->modDownTo(cipher, r.logq);
        scheme->multAndEqual(w, r);
        scheme->reScaleByAndEqual(w, logp);

        Ciphertext rSquare = scheme->square(r);
        scheme->reScaleByAndEqual(rSquare, logp);

        scheme->multAndEqual(w, rSquare);
        scheme->reScaleByAndEqual(w, logp);

        nprint("w" + to_string(i), w);

        // r <- (3r - w) / 2
        scheme->multByConstAndEqual(r, 3.0, 0);
        scheme->modDownToAndEqual(r, w.logq);
        scheme->negateAndEqual(w);
        scheme->addAndEqual(r, w);
        scheme->divByPo2AndEqual(r, 1);

        nprint("r" + to_string(i), r);
        // cout << "r.logq = " << r.logq << endl;
    }

    scheme->checkLevelAndBoot(r, 1, *bootHelper, param);
    scheme->checkAndBoot(cipher, cipher.logq < r.logq, *bootHelper, param);

    scheme->modDownToAndEqualModified(cipher, r, *bootHelper, param); // a - 1
    scheme->multAndEqualWithBoot(cipher, r, *bootHelper, param);
    scheme->reScaleByAndEqual(cipher, logp);
}

void EncAlgo::minMax(Ciphertext& minCipher, Ciphertext& maxCipher) {
    Ciphertext x = scheme->add(minCipher, maxCipher);
    Ciphertext y = scheme->sub(minCipher, maxCipher);

    scheme->divByPo2AndEqual(x, 1); // x - logp + 1

    scheme->checkModulusAndBoot(y, param.logp + 1, *bootHelper, param); //FIND_HERE
    Ciphertext yBefore = y;    

    nprint("before square, y", y);
    scheme->resetImagErrorAndEqual(y);
    scheme->squareAndEqual(y);
    scheme->reScaleByAndEqual(y, param.logp);
    scheme->resetImagErrorAndEqual(y);   
    nprint("after square, y", y);

    sqrtAlgorithm(y);

    if(printCondition) {
        complex<double>* yBefDec = scheme->decrypt(*secretKey, yBefore);    
        complex<double>* yAftDec = scheme->decrypt(*secretKey, y);
        complex<double>* yBefAbs = new complex<double>[y.n];
        for (int i = 0; i < y.n; i++) {
            yBefAbs[i] = abs(yBefDec[i]);
        }
        cout << "after approxSqrt" << endl;
        cout << "sqrt(y) (logQ = " << yBefore.logq << ") // approxsqrt(y) (logQ = " << y.logq << ")" << endl;
        PrintUtils::printFewArrays(yBefDec, yAftDec, y.n);
        PrintUtils::averageDifference(yBefAbs, yAftDec, y.n);
    }
    
    scheme->divByPo2AndEqual(y, 1);

    scheme->modDownToAndEqualModified(x, y, *bootHelper, param);

    maxCipher = scheme->add(x, y);
    minCipher = scheme->sub(x, y);
}

void EncAlgo::newMinMax(Ciphertext& minCipher, Ciphertext& maxCipher) {
    // Ciphertext max = scheme->sub(maxCipher, minCipher);
    // evalFcn(max);
    // scheme->addConstAndEqual(max, 1.0, param.logp);
    // scheme->divByPo2AndEqual(max, 1);
    // Ciphertext min = scheme->negate(max);
    // scheme->addConstAndEqual(min, 1.0, param.logp);
    
    // scheme->modDownToAndEqualModified(minCipher, max, *bootHelper, param);
    // scheme->modDownToAndEqualModified(maxCipher, max, *bootHelper, param);
    // Ciphertext minMin = scheme->mult(minCipher, min);
    // Ciphertext minMax = scheme->mult(minCipher, max);
    // Ciphertext maxMin = scheme->mult(maxCipher, min);
    // Ciphertext maxMax = scheme->mult(maxCipher, max);

    // minCipher = scheme->add(minMax, maxMin);
    // maxCipher = scheme->add(minMin, maxMax);
    // scheme->reScaleByAndEqual(minCipher, param.logp);
    // scheme->reScaleByAndEqual(maxCipher, param.logp);

    Ciphertext add = scheme->add(maxCipher, minCipher);
    Ciphertext sub = scheme->sub(maxCipher, minCipher);
    Ciphertext minMax = sub;
    
    nprint("a - b", minMax);
    evalFcn(minMax);
    nprint("f(a - b)", minMax);

    scheme->divByPo2AndEqual(add, 1);
    scheme->divByPo2AndEqual(sub, 1);
    scheme->modDownToAndEqualModified(minMax, sub, *bootHelper, param);
    scheme->multAndEqual(sub, minMax);
    scheme->reScaleByAndEqual(sub, param.logp);

    nprint("0 < (a-b)/ * f(a-b)", sub);

    scheme->modDownToAndEqualModified(add, sub, *bootHelper, param);
    
    maxCipher = scheme->add(add, sub);
    minCipher = scheme->sub(add, sub);
}

void EncAlgo::encSwap(Ciphertext& cipher, ZZ* mask, long dist, bool increase) {
    long n = cipher.n;
    // scheme->checkModulusAndBoot(cipher, param.logp + 1, *bootHelper, param);
    scheme->checkModulusAndBoot(cipher, 3 * param.logp, *bootHelper, param);

    nprint("before reset imag, cipher", cipher);
    
    
    Ciphertext dummy = cipher;
    scheme->multByPolyAndEqualWithBoot(dummy, mask, *bootHelper, param);
    scheme->reScaleByAndEqual(dummy, param.logp);
    scheme->modDownToAndEqualModified(cipher, dummy, *bootHelper, param);
    scheme->subAndEqual(cipher, dummy);

    scheme->addByPolyAndEqual(cipher, mask, param.logp);
    scheme->rightRotateAndEqualConditional(dummy, dist, increase);

    nprint("before minmax, cipher", cipher);
    nprint("before minmax, dummy", dummy);

    minMaxAlgorithm(dummy, cipher);

    nprint("after minmax, cipher", cipher);
    nprint("after minmax, dummy", dummy);

    scheme->subByPolyAndEqual(cipher, mask, param.logp);
    scheme->leftRotateAndEqualConditional(dummy, dist, increase); 
    
    scheme->addAndEqual(cipher, dummy);   
}

void EncAlgo::selfBitonicMerge(Ciphertext& cipher, ZZ** mask, bool increase) {
    for(int i = 0; i < param.log2n; i++) {
        long logDist = param.log2n - 1 - i;
        nprint("                   - selfBitonicMerge, " + to_string(logDist));
        encSwap(cipher, mask[i], 1 << logDist, increase);
    }
}

void EncAlgo::selfTableMerge(Ciphertext& cipher, ZZ** mask, bool increase) {
    for(int i = 0; i < param.log2n - logDataNum; i++) {
        nprint("                   - selfBitonicMerge, " + to_string(param.log2n - 1 - i));
        // encSwapTable(cipher, mask[i], mask[i], maskMergeTablePoly[inc][i], maskMergeTableColPoly[inc][i], dist, increase);
        // compAndSwapTable(cipher, logDataNum, colNum, mask[0][i], mask[1][i], mask[2][i], mask[3][i], 1 << (param.log2n - 1 - i), scheme, ring, bootHelper, sk);
    }
}

void EncAlgo::reverse(Ciphertext& cipher, ZZ** mask) {
    scheme->leftRotateFastAndEqual(cipher, 1 << (param.log2n - 1));
    for(int i = 1; i < param.log2n; i++) {
        Ciphertext dummy = cipher;
        scheme->multByPolyAndEqualWithBoot(dummy, mask[i], *bootHelper, param);
        scheme->reScaleByAndEqual(dummy, param.logp);
        scheme->modDownToAndEqualModified(cipher, dummy, *bootHelper, param);
        scheme->subAndEqual(cipher, dummy);

        scheme->leftRotateFastAndEqual(cipher, 1 << (param.log2n - i - 1));
        scheme->rightRotateFastAndEqual(dummy, 1 << (param.log2n - i - 1));
        scheme->addAndEqual(cipher, dummy);
    }
}

void EncAlgo::reverse(Ciphertext& cipher, ZZ** maskLeft, ZZ** maskRight, long level, bool increase) {
    for(int i = 1; i <= level; i++) {
        nprint("i = " + to_string(i));
        scheme->checkModulusAndBoot(cipher, param.logp + 1, *bootHelper, param);
        scheme->resetImagErrorAndEqual(cipher);

        Ciphertext left = scheme->multByPoly(cipher, maskLeft[i], param.logp);
        Ciphertext right = scheme->multByPoly(cipher, maskRight[i], param.logp);
        scheme->reScaleByAndEqual(left, param.logp);
        scheme->reScaleByAndEqual(right, param.logp);
        scheme->modDownToAndEqualModified(cipher, left, *bootHelper, param);
        scheme->subAndEqual(cipher, left);
        scheme->subAndEqual(cipher, right);
        nprint("rotate by " + to_string(level - i));
        scheme->rightRotateAndEqualConditional(left, 1 << (level - i), increase);
        scheme->leftRotateAndEqualConditional(right, 1 << (level - i), increase);
        scheme->addAndEqual(cipher, left);
        scheme->addAndEqual(cipher, right);
    }
}

void EncAlgo::halfCleaner(Ciphertext& cipher, ZZ* mask, long dist, bool increase) {
    scheme->checkLevelAndBoot(cipher, 1, *bootHelper, param);
        
    Ciphertext dummy = scheme->multByPoly(cipher, mask, param.logp);
    scheme->reScaleByAndEqual(dummy, param.logp);
    scheme->modDownByAndEqual(cipher, param.logp);
    scheme->subAndEqual(cipher, dummy);
    

    scheme->rightRotateAndEqualConditional(dummy, dist, increase);
    scheme->addByPolyAndEqual(cipher, mask, param.logp);

    minMaxAlgorithm(dummy, cipher);

    scheme->leftRotateAndEqualConditional(dummy, dist, increase);
    scheme->subByPolyAndEqual(cipher, mask, param.logp);
    
    scheme->addAndEqual(cipher, dummy);
}

void EncAlgo::approxInverse(Ciphertext& cipher) {
    scheme->negateAndEqual(cipher);

    Ciphertext a = scheme->addConst(cipher, 2.0, param.logp);
    Ciphertext b = scheme->addConst(cipher, 1.0, param.logp);

    for (int _ = 0; _ < iter[invLoc]; _++) {
        scheme->squareAndEuqalWithBoot(b, *bootHelper, param);
        scheme->reScaleByAndEqual(b, param.logp);
        Ciphertext bPlusOne = scheme->addConst(b, 1.0, param.logp);
        scheme->modDownToAndEqualModified(a, bPlusOne, *bootHelper, param);
        nprint("start a", a);
        nprint("bPlusOne", bPlusOne);
        
        // scheme->multAndEqualWithBoot(a, bPlusOne, *bootHelper, param);
        if (a.logq - param.logp < param.logq) {
            bootHelper->bootstrapping(a, param.logq, param.logQ, param.logT);
            bootHelper->bootstrapping(bPlusOne, param.logq, param.logQ, param.logT);
            nprint("after boot a", a);
            nprint("after boot bPlusOne", bPlusOne);
        }
        
        scheme->multAndEqual(a, bPlusOne);
        scheme->reScaleByAndEqual(a, param.logp);
        nprint("mult a", a);
        scheme->modDownToAndEqualModified(b, a, *bootHelper, param);
    }
    cipher = a;
    nprint("last a", cipher);
}

void EncAlgo::comparison(Ciphertext& a, Ciphertext& b) {
    /*
        uses approxInv * (compIter + 1) times
    */
    // Input messages are in (0, 1)
    scheme->addConstAndEqual(a, 0.5, param.logp);
    scheme->addConstAndEqual(b, 0.5, param.logp); // now in (1/2, 3/2)

    // normalize cipherLeft, cipher2
    scheme->checkAndBoot(a, a.logq - 1 < param.logq, *bootHelper, param);
    scheme->checkAndBoot(b, b.logq < a.logq, *bootHelper, param);
    
    // sum <- 2/(a+b) = Inv((a+b)/2)
    Ciphertext sum = scheme->add(a, b);
    scheme->divByPo2AndEqual(sum, 1);

    nprint("bofore inverse : sum", sum);

    approxInverse(sum);

    nprint("after inverse : sum", sum);

    
    // cout << "a = " << a.logq << ", inv(sum) = " << sum.logq << endl;

    // to compute a <- (a/2) * sum = a/(a+b)
    scheme->checkAndBoot(a, a.logq - param.logp - 1 < param.logq, *bootHelper, param);
    scheme->checkAndBoot(sum, sum.logq - param.logp - 1 < param.logq, *bootHelper, param);
    
    scheme->divByPo2AndEqual(a, 1);
    scheme->modDownToAndEqualModified(a, sum, *bootHelper, param);
    scheme->multAndEqualWithBoot(a, sum, *bootHelper, param);
    scheme->reScaleByAndEqual(a, param.logp);

    nprint("mult by sum, a", a);


    for (int _ = 0; _ < iter[compLoc]; _++) {
        // check level for a^m
        scheme->checkLevelAndBoot(a, 2, *bootHelper, param);

        // update b <- 1 - a
        b = scheme->negate(a);
        scheme->addConstAndEqual(b, 1.0, param.logp);

        nprint("before squaring, a", a);
        nprint("before squaring, b", b);

        // Fix m = 4
        for (int i = 0; i < 2; i++) { // consumes 2 levels
            scheme->squareAndEuqalWithBoot(a, *bootHelper, param);
            scheme->reScaleByAndEqual(a, param.logp);
            
            scheme->squareAndEuqalWithBoot(b, *bootHelper, param);
            scheme->reScaleByAndEqual(b, param.logp);
        }
        nprint("after squaring, a", a);
        nprint("after squaring, b", b);
        
        Ciphertext inv = scheme->add(a, b);

        nprint("before inv, inv", inv);
        nprint("inv.logq = " + to_string(inv.logq));
        if(inv.logq - (iter[invLoc] + 2) * param.logp < param.logq) {
            scheme->checkAndBoot(inv, true, *bootHelper,param);
            nprint("Bootstrap inv, inv", inv);
            nprint("inv.logq = " + to_string(inv.logq));
        }
        approxInverse(inv); // consumses invIter + 1 levels

        nprint("after inv, inv", inv);
        nprint("inv.logq = " + to_string(inv.logq));

        nprint("a.logq = " + to_string(a.logq) + ", inv.logq = " + to_string(inv.logq));
        
        if(a.logq - param.logp < param.logq || inv.logq - param.logp < param.logq) {
            scheme->checkAndBoot(a, true, *bootHelper, param);
            scheme->checkAndBoot(inv, true, *bootHelper, param);
        }

        scheme->modDownToAndEqualModified(a, inv, *bootHelper, param);

        nprint("a.logq = " + to_string(a.logq) + ", inv.logq = " + to_string(inv.logq));

        scheme->multAndEqualWithBoot(a, inv, *bootHelper, param);
        scheme->reScaleByAndEqual(a, param.logp);

        nprint("Iteration : final a", a);
    }
    b = scheme->negate(a);
    scheme->addConstAndEqual(b, 1.0, param.logp);
}

void EncAlgo::newComparison(Ciphertext& a, Ciphertext& b) {
    Ciphertext max = scheme->sub(a, b);
    evalFcn(max);
    scheme->addConstAndEqual(max, 1.0, param.logp);
    scheme->divByPo2AndEqual(max, 1);
    Ciphertext min = scheme->negate(max);
    scheme->addConstAndEqual(min, 1.0, param.logp);

    a = max;
    b = min;
}

void EncAlgo::minMaxTable(Ciphertext& minCipher, Ciphertext& maxCipher, Ciphertext& minCipherTable, Ciphertext& maxCipherTable, ZZ* mask, ZZ* maskTable) {
    /*
        Inputs:
            minCipher : (a, a1, a2, a3), ( ... )
            maxCipher : (b, b1, b2, b3), ( ... )
            minCipherTable : (a, 0, 0, 0), ( ... )
            maxCipherTable : (b, 0, 0, 0), ( ... )
    */
    comparisonAlgorithm(minCipherTable, maxCipherTable);
    /*
        after comp :
            minCipherTable : (comp(a,b), *, *, *), ( ... )
            maxCipherTable : (comp(b,a), *, *, *), ( ... )
    */

    nprint("after comp");
    nprint("minTable", minCipherTable);
    nprint("maxTable", maxCipherTable);

    scheme->checkLevelAndBoot(minCipherTable, 2, *bootHelper, param);
    scheme->checkLevelAndBoot(maxCipherTable, 2, *bootHelper, param);
    scheme->checkLevelAndBoot(minCipher, 1, *bootHelper, param);
    scheme->checkLevelAndBoot(maxCipher, 1, *bootHelper, param);

    // mult maskTablePoly = (1, 0, 0, 0) to remove * parts
    scheme->multByPolyAndEqualWithBoot(minCipherTable, maskTable, *bootHelper, param);
    scheme->multByPolyAndEqualWithBoot(maxCipherTable, maskTable, *bootHelper, param);
    scheme->reScaleByAndEqual(minCipherTable, param.logp);
    scheme->reScaleByAndEqual(maxCipherTable, param.logp);

    /*
        after comp :
            minCipherTable : (comp(a,b), 0, 0, 0), ( ... )
            maxCipherTable : (comp(b,a), 0, 0, 0), ( ... )
    */

    nprint("after masking");
    nprint("minTable", minCipherTable);
    nprint("maxTable", maxCipherTable);


    /*
        copy :
            minCipherTable : (comp(a,b), comp(a,b), comp(a,b), comp(a,b)), ( ... )
            maxCipherTable : (comp(b,a), comp(b,a), comp(b,a), comp(b,a)), ( ... )
    */
    scheme->leftRotateFastAndEqual(minCipherTable, colNum);
    scheme->leftRotateFastAndEqual(maxCipherTable, colNum);

    
    for (int i = 0; i < logDataNum; i++) {
        Ciphertext tmpMinTable = scheme->rightRotateFast(minCipherTable, 1 << i);
        Ciphertext tmpMaxTable = scheme->rightRotateFast(maxCipherTable, 1 << i);
        scheme->addAndEqual(minCipherTable, tmpMinTable);
        scheme->addAndEqual(maxCipherTable, tmpMaxTable);
    }

    nprint("after rot");
    nprint("minTable", minCipherTable);
    nprint("maxTable", maxCipherTable);

    // Flip <- 1 - Table
    Ciphertext minCipherTableFlip = scheme->negate(minCipherTable);
    Ciphertext maxCipherTableFlip = scheme->negate(maxCipherTable);

    scheme->addByPolyAndEqual(minCipherTableFlip, mask, param.logp);
    scheme->addByPolyAndEqual(maxCipherTableFlip, mask, param.logp);
    
    // scheme->addByPolyAndEqual(cipherTableFlip, maskPoly, param.logp);
    // scheme->addByPolyAndEqual(cipherTableFlipRight, maskRightPoly, param.logp);


    // scheme->decryptAndPrint("tableFlip", secretKey, minCipherTableFlip);
    // scheme->decryptAndPrint("tableFlipRight", secretKey, maxCipherTableFlip);

    if (minCipher.logq < minCipherTable.logq) {
        scheme->modDownToAndEqual(minCipherTable, minCipher.logq);
        scheme->modDownToAndEqual(minCipherTableFlip, minCipher.logq);
        scheme->modDownToAndEqual(maxCipherTable, minCipher.logq);
        scheme->modDownToAndEqual(maxCipherTableFlip, minCipher.logq);
    } else {
        scheme->modDownToAndEqual(minCipher, minCipherTable.logq);
        scheme->modDownToAndEqual(maxCipher, minCipherTable.logq);
    }
        
    Ciphertext minCipherSmall = scheme->multWithBoot(minCipher, minCipherTableFlip, *bootHelper, param);
    Ciphertext minCipherBig = scheme->multWithBoot(minCipher, minCipherTable, *bootHelper, param);
    Ciphertext maxCipherSmall = scheme->multWithBoot(maxCipher, maxCipherTableFlip, *bootHelper, param);
    Ciphertext maxCipherBig = scheme->multWithBoot(maxCipher, maxCipherTable, *bootHelper, param);
    scheme->reScaleByAndEqual(minCipherSmall, param.logp);
    scheme->reScaleByAndEqual(minCipherBig, param.logp);
    scheme->reScaleByAndEqual(maxCipherSmall, param.logp);
    scheme->reScaleByAndEqual(maxCipherBig, param.logp);

    nprint("check");
    nprint("minSmall", minCipherSmall);
    nprint("minBig", minCipherBig);
    nprint("maxSmall", maxCipherSmall);
    nprint("maxBig", maxCipherBig);

    minCipher = scheme->add(minCipherSmall, maxCipherSmall);
    maxCipher = scheme->add(minCipherBig, maxCipherBig);

    nprint("result");
    nprint("minCipher", minCipher);
    nprint("maxCipher", maxCipher);
}

void EncAlgo::encSwapTable(Ciphertext& cipher, ZZ* maskLeft, ZZ* maskRight, ZZ* maskTableLeft, ZZ* maskTableRight, long dist, bool increase) {
    // cout << "start compAndSwapTable with dist = " << dist << endl; 
    scheme->checkLevelAndBoot(cipher, 3, *bootHelper, param);

    Ciphertext cipherLeft = scheme->multByPolyWithBoot(cipher, maskLeft, *bootHelper, param);
    Ciphertext cipherRight = scheme->multByPolyWithBoot(cipher, maskRight, *bootHelper, param);
    Ciphertext cipherTable = scheme->multByPolyWithBoot(cipher, maskTableLeft, *bootHelper, param);
    Ciphertext cipherTableRight = scheme->multByPolyWithBoot(cipher, maskTableRight, *bootHelper, param); 
    scheme->reScaleByAndEqual(cipherLeft, param.logp);
    scheme->reScaleByAndEqual(cipherRight, param.logp);
    scheme->reScaleByAndEqual(cipherTable, param.logp);
    scheme->reScaleByAndEqual(cipherTableRight, param.logp);

    nprint("cipherLeft", cipherLeft);
    nprint("cipherRight", cipherRight);
    nprint("cipherTableLeft", cipherTable);
    nprint("cipherTableRight", cipherTableRight);

    scheme->modDownToAndEqual(cipher, cipherLeft.logq);
    scheme->subAndEqual(cipher, cipherLeft);
    scheme->subAndEqual(cipher, cipherRight);

    scheme->rightRotateAndEqualConditional(cipherLeft, dist, increase);
    scheme->rightRotateAndEqualConditional(cipherTable, dist, increase);

    minMaxTable(cipherLeft, cipherRight, cipherTable, cipherTableRight, maskRight, maskTableRight);

    nprint("\n.......after minMax.........\n");
    nprint("cipherLeft", cipherLeft);
    nprint("cipherRight", cipherRight);
    nprint("cipherTableLeft", cipherTable);
    nprint("cipherTableRight", cipherTableRight);


    scheme->leftRotateAndEqualConditional(cipherLeft, dist, increase);

    scheme->modDownToAndEqual(cipher, cipherLeft.logq);

    scheme->addAndEqual(cipher, cipherLeft);
    scheme->addAndEqual(cipher, cipherRight);
}