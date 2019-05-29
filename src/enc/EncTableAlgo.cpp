// #include "EncTableAlgo.h"

// EncTableAlgo::EncTableAlgo(Parameter _param, long _logDataNum, long _colNum, long _invIter, long _compIter, bool _printCondition) :
//         EncInterface(_param, _printCondition),
//         logDataNum(_logDataNum), colNum(_colNum),
//         invIter(_invIter), compIter(_compIter) {
//     std::cout << "***************************" << std::endl;
//     std::cout << "Other Parameters" << std::endl;
//     std::cout << "logDataNum = " << logDataNum << ", colNum = " << colNum << std::endl;
//     std::cout << "invIter = " << invIter << ", compIter = " << compIter << std::endl;
//     std::cout << "want to print? : " << printCondition << std::endl;
//     std::cout << "***************************" << std::endl;
//     std::cout << std::endl;

// }

// void EncTableAlgo::approxInverse(Ciphertext& cipher) {
//     // cout << "Asdf " << endl;

//     // scheme->negateAndEqual(cipher);
//     // cout << "Asdf " << endl;

//     // Ciphertext a = scheme->addConst(cipher, 2.0, param.logp);
//     // Ciphertext b = scheme->addConst(cipher, 1.0, param.logp);
//     // cout << "Asdf " << endl;
    
//     for (int _ = 0; _ < iter[invLoc]; _++) {
//         scheme->squareAndEuqalWithBoot(b, *bootHelper, param);
//         scheme->reScaleByAndEqual(b, param.logp);
//         Ciphertext bPlusOne = scheme->addConst(b, 1.0, param.logp);
//         scheme->modDownToAndEqualModified(a, bPlusOne, *bootHelper, param);
//         nprint("start a", a);
//         nprint("bPlusOne", bPlusOne);
        
//         // scheme->multAndEqualWithBoot(a, bPlusOne, *bootHelper, param);
//         if (a.logq - param.logp < param.logq) {
//             bootHelper->bootstrapping(a, param.logq, param.logQ, param.logT);
//             bootHelper->bootstrapping(bPlusOne, param.logq, param.logQ, param.logT);
//             nprint("after boot a", a);
//             nprint("after boot bPlusOne", bPlusOne);
//         }
        
//         scheme->multAndEqual(a, bPlusOne);
//         scheme->reScaleByAndEqual(a, param.logp);
//         nprint("mult a", a);
//         scheme->modDownToAndEqualModified(b, a, *bootHelper, param);
//     }
//     cipher = a;
//     nprint("last a", cipher);
// }

// void EncTableAlgo::comparison(Ciphertext& a, Ciphertext& b) {
//     /*
//         uses approxInv * (compIter + 1) times
//     */
//     // Input messages are in (0, 1)
//     scheme->addConstAndEqual(a, 0.5, param.logp);
//     scheme->addConstAndEqual(b, 0.5, param.logp); // now in (1/2, 3/2)

//     // normalize cipherLeft, cipher2
//     scheme->checkAndBoot(a, a.logq - 1 < param.logq, *bootHelper, param);
//     scheme->checkAndBoot(b, b.logq < a.logq, *bootHelper, param);
    
//     // sum <- 2/(a+b) = Inv((a+b)/2)
//     Ciphertext sum = scheme->add(a, b);
//     scheme->divByPo2AndEqual(sum, 1);

//     nprint("bofore inverse : sum", sum);

//     approxInverse(sum);

//     nprint("after inverse : sum", sum);

    
//     // cout << "a = " << a.logq << ", inv(sum) = " << sum.logq << endl;

//     // to compute a <- (a/2) * sum = a/(a+b)
//     scheme->checkAndBoot(a, a.logq - param.logp - 1 < param.logq, *bootHelper, param);
//     scheme->checkAndBoot(sum, sum.logq - param.logp - 1 < param.logq, *bootHelper, param);
    
//     scheme->divByPo2AndEqual(a, 1);
//     scheme->modDownToAndEqualModified(a, sum, *bootHelper, param);
//     scheme->multAndEqualWithBoot(a, sum, *bootHelper, param);
//     scheme->reScaleByAndEqual(a, param.logp);

//     nprint("mult by sum, a", a);


//     for (int _ = 0; _ < compIter; _++) {
//         // check level for a^m
//         scheme->checkLevelAndBoot(a, 2, *bootHelper, param);

//         // update b <- 1 - a
//         b = scheme->negate(a);
//         scheme->addConstAndEqual(b, 1.0, param.logp);

//         nprint("before squaring, a", a);
//         nprint("before squaring, b", b);

//         // Fix m = 4
//         for (int i = 0; i < 2; i++) { // consumes 2 levels
//             scheme->squareAndEuqalWithBoot(a, *bootHelper, param);
//             scheme->reScaleByAndEqual(a, param.logp);
            
//             scheme->squareAndEuqalWithBoot(b, *bootHelper, param);
//             scheme->reScaleByAndEqual(b, param.logp);
//         }
//         nprint("after squaring, a", a);
//         nprint("after squaring, b", b);
        
//         Ciphertext inv = scheme->add(a, b);

//         nprint("before inv, inv", inv);
//         cout << "inv.logq = " << inv.logq << endl;
//         if(inv.logq - (invIter + 2) * param.logp < param.logq) {
//             scheme->checkAndBoot(inv, true, *bootHelper,param);
//             nprint("Bootstrap inv, inv", inv);
//             cout << "inv.logq = " << inv.logq << endl;
//         }
//         approxInverse(inv); // consumses invIter + 1 levels

//         nprint("after inv, inv", inv);
//         cout << "inv.logq = " << inv.logq << endl;

//         cout << "a = " << a.logq << ", inv = " << inv.logq << endl;    

        
//         if(a.logq - param.logp < param.logq || inv.logq - param.logp < param.logq) {
//             scheme->checkAndBoot(a, true, *bootHelper, param);
//             scheme->checkAndBoot(inv, true, *bootHelper, param);
//         }

//         scheme->modDownToAndEqualModified(a, inv, *bootHelper, param);

//         cout << "a = " << a.logq << ", inv = " << inv.logq << endl;    

//         scheme->multAndEqualWithBoot(a, inv, *bootHelper, param);
//         scheme->reScaleByAndEqual(a, param.logp);

//         nprint("Iteration : final a", a);
//     }
//     b = scheme->negate(a);
//     scheme->addConstAndEqual(b, 1.0, param.logp);
//     // cout << a.logq << ", " << b.logq << endl;    
//     // cout << "** End Comparison ** " << endl;
// }

// void EncTableAlgo::minMaxTable(Ciphertext& minCipher, Ciphertext& maxCipher, Ciphertext& minCipherTable, Ciphertext& maxCipherTable, ZZ* mask, ZZ* maskTable) {
//     /*
//         Inputs:
//             minCipher : (a, a1, a2, a3), ( ... )
//             maxCipher : (b, b1, b2, b3), ( ... )
//             minCipherTable : (a, 0, 0, 0), ( ... )
//             maxCipherTable : (b, 0, 0, 0), ( ... )
//     */
//     comparison(minCipherTable, maxCipherTable);
//     /*
//         after comp :
//             minCipherTable : (comp(a,b), *, *, *), ( ... )
//             maxCipherTable : (comp(b,a), *, *, *), ( ... )
//     */

//     nprint("after comp");
//     nprint("minTable", minCipherTable);
//     nprint("maxTable", maxCipherTable);

//     scheme->checkLevelAndBoot(minCipherTable, 2, *bootHelper, param);
//     scheme->checkLevelAndBoot(maxCipherTable, 2, *bootHelper, param);
//     scheme->checkLevelAndBoot(minCipher, 1, *bootHelper, param);
//     scheme->checkLevelAndBoot(maxCipher, 1, *bootHelper, param);

//     // mult maskTablePoly = (1, 0, 0, 0) to remove * parts
//     scheme->multByPolyAndEqualWithBoot(minCipherTable, maskTable, *bootHelper, param);
//     scheme->multByPolyAndEqualWithBoot(maxCipherTable, maskTable, *bootHelper, param);
//     scheme->reScaleByAndEqual(minCipherTable, param.logp);
//     scheme->reScaleByAndEqual(maxCipherTable, param.logp);

//     /*
//         after comp :
//             minCipherTable : (comp(a,b), 0, 0, 0), ( ... )
//             maxCipherTable : (comp(b,a), 0, 0, 0), ( ... )
//     */

//     nprint("after masking");
//     nprint("minTable", minCipherTable);
//     nprint("maxTable", maxCipherTable);


//     /*
//         copy :
//             minCipherTable : (comp(a,b), comp(a,b), comp(a,b), comp(a,b)), ( ... )
//             maxCipherTable : (comp(b,a), comp(b,a), comp(b,a), comp(b,a)), ( ... )
//     */
//     scheme->leftRotateFastAndEqual(minCipherTable, colNum);
//     scheme->leftRotateFastAndEqual(maxCipherTable, colNum);

    
//     for (int i = 0; i < logDataNum; i++) {
//         Ciphertext tmpMinTable = scheme->rightRotateFast(minCipherTable, 1 << i);
//         Ciphertext tmpMaxTable = scheme->rightRotateFast(maxCipherTable, 1 << i);
//         scheme->addAndEqual(minCipherTable, tmpMinTable);
//         scheme->addAndEqual(maxCipherTable, tmpMaxTable);
//     }

//     nprint("after rot");
//     nprint("minTable", minCipherTable);
//     nprint("maxTable", maxCipherTable);

//     // Flip <- 1 - Table
//     Ciphertext minCipherTableFlip = scheme->negate(minCipherTable);
//     Ciphertext maxCipherTableFlip = scheme->negate(maxCipherTable);

//     scheme->addByPolyAndEqual(minCipherTableFlip, mask, param.logp);
//     scheme->addByPolyAndEqual(maxCipherTableFlip, mask, param.logp);
    
//     // scheme->addByPolyAndEqual(cipherTableFlip, maskPoly, param.logp);
//     // scheme->addByPolyAndEqual(cipherTableFlipRight, maskRightPoly, param.logp);


//     // scheme->decryptAndPrint("tableFlip", secretKey, minCipherTableFlip);
//     // scheme->decryptAndPrint("tableFlipRight", secretKey, maxCipherTableFlip);

//     if (minCipher.logq < minCipherTable.logq) {
//         scheme->modDownToAndEqual(minCipherTable, minCipher.logq);
//         scheme->modDownToAndEqual(minCipherTableFlip, minCipher.logq);
//         scheme->modDownToAndEqual(maxCipherTable, minCipher.logq);
//         scheme->modDownToAndEqual(maxCipherTableFlip, minCipher.logq);
//     } else {
//         scheme->modDownToAndEqual(minCipher, minCipherTable.logq);
//         scheme->modDownToAndEqual(maxCipher, minCipherTable.logq);
//     }
        
//     Ciphertext minCipherSmall = scheme->multWithBoot(minCipher, minCipherTableFlip, *bootHelper, param);
//     Ciphertext minCipherBig = scheme->multWithBoot(minCipher, minCipherTable, *bootHelper, param);
//     Ciphertext maxCipherSmall = scheme->multWithBoot(maxCipher, maxCipherTableFlip, *bootHelper, param);
//     Ciphertext maxCipherBig = scheme->multWithBoot(maxCipher, maxCipherTable, *bootHelper, param);
//     scheme->reScaleByAndEqual(minCipherSmall, param.logp);
//     scheme->reScaleByAndEqual(minCipherBig, param.logp);
//     scheme->reScaleByAndEqual(maxCipherSmall, param.logp);
//     scheme->reScaleByAndEqual(maxCipherBig, param.logp);

//     nprint("check");
//     nprint("minSmall", minCipherSmall);
//     nprint("minBig", minCipherBig);
//     nprint("maxSmall", maxCipherSmall);
//     nprint("maxBig", maxCipherBig);

//     minCipher = scheme->add(minCipherSmall, maxCipherSmall);
//     maxCipher = scheme->add(minCipherBig, maxCipherBig);

//     nprint("result");
//     nprint("minCipher", minCipher);
//     nprint("maxCipher", maxCipher);
// }

// void EncTableAlgo::compAndSwapTable(Ciphertext& cipher, long logDataNum, long colNum, ZZ* maskLeft, ZZ* maskRight, ZZ* maskTableLeft, ZZ* maskTableRight, long dist, bool increase) {
//     // cout << "start compAndSwapTable with dist = " << dist << endl; 
//     scheme->checkLevelAndBoot(cipher, 3, *bootHelper, param);

//     Ciphertext cipherLeft = scheme->multByPolyWithBoot(cipher, maskLeft, *bootHelper, param);
//     Ciphertext cipherRight = scheme->multByPolyWithBoot(cipher, maskRight, *bootHelper, param);
//     Ciphertext cipherTable = scheme->multByPolyWithBoot(cipher, maskTableLeft, *bootHelper, param);
//     Ciphertext cipherTableRight = scheme->multByPolyWithBoot(cipher, maskTableRight, *bootHelper, param); 
//     scheme->reScaleByAndEqual(cipherLeft, param.logp);
//     scheme->reScaleByAndEqual(cipherRight, param.logp);
//     scheme->reScaleByAndEqual(cipherTable, param.logp);
//     scheme->reScaleByAndEqual(cipherTableRight, param.logp);

//     scheme->modDownToAndEqual(cipher, cipherLeft.logq);
//     scheme->subAndEqual(cipher, cipherLeft);
//     scheme->subAndEqual(cipher, cipherRight);

//     scheme->rightRotateAndEqualConditional(cipherLeft, dist, increase);
//     scheme->rightRotateAndEqualConditional(cipherTable, dist, increase);

//     minMaxTable(cipherLeft, cipherRight, cipherTable, cipherTableRight, maskLeft, maskTableLeft);    

//     scheme->leftRotateAndEqualConditional(cipherLeft, dist, increase);

//     scheme->modDownToAndEqual(cipher, cipherLeft.logq);

//     scheme->addAndEqual(cipher, cipherLeft);
//     scheme->addAndEqual(cipher, cipherRight);
// }