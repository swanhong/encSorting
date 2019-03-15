#include "../HEAAN/src/HEAAN.h"
#include "../bootsrc/new_bootstrapping.h"
#include "Parameter.h"

void multWithBootAndEqual(Scheme& scheme, Ciphertext& cipher1, Ciphertext& cipher2, BootHelper boothelper, Parameter parameter);

void fcnEncSqrtWithBoot(Ciphertext& outCipher, const Ciphertext& inCipher, Parameter parameter, long d, Scheme& scheme, BootHelper& boothelper);

void fcnEncMaxMinWithBoot(Ciphertext& maxCipher, Ciphertext& minCipher, Ciphertext& input1, Ciphertext& input2, Parameter parameter, long iter, Scheme& scheme, BootHelper& boothelper);

void fcnEncCompAndSwapWithBoot(Ciphertext& outCipher, const Ciphertext& inCipher, double* mask, long dist, Parameter parameter, long iter, Scheme& scheme, BootHelper& boothelper);

long encBatcherOddEvenSortRec(Ciphertext& sortedCipher, const Ciphertext& inCipher, long logNum, long logJump, long loc, Parameter parameter, long iter, double** mask, Scheme& scheme, BootHelper& boothelper);

void encBatcherOddEvenSort(Ciphertext& sortedCipher, const Ciphertext& inCipher, Parameter parameter, long iter, Scheme& scheme, BootHelper& boothelper);

void encBatcherOddEvenSortWithDecrypt(Ciphertext& sortedCipher, const Ciphertext& inCipher, Parameter parameter, long iter, Scheme& scheme, BootHelper& boothelper, SecretKey& secretKey);

long encBatcherOddEvenSortRecWithDecrypt(Ciphertext& sortedCipher, const Ciphertext& inCipher, long logNum, long logJump, long loc, Parameter parameter, long iter, double** mask, Scheme& scheme, BootHelper& boothelper, SecretKey& secretKey);

void fcnEncCompAndSwapWithBootAndDecrypt(Ciphertext& outCipher, Ciphertext& inCipher, double* mask, long dist, Parameter parameter, long iter, Scheme& scheme, BootHelper& boothelper, SecretKey& secretKey);

void fcnEncMaxMinWithBootAndDecrypt(Ciphertext& maxCipher, Ciphertext& minCipher, Ciphertext& input1, Ciphertext& input2, Parameter parameter, long iter, Scheme& scheme, BootHelper& boothelper, SecretKey& secretKey);

void fcnEncSqrtWithBootWithDecrypt(Ciphertext& outCipher, Ciphertext& inCipher, Parameter parameter, long d, Scheme& scheme, BootHelper& boothelper, SecretKey& secretKey);
// ======================================================================
// ======================================================================
// ======================================================================
// ======================================================================
// ======================================================================

void multWithBootAndEqual(Scheme& scheme, Ciphertext& cipher1, Ciphertext& cipher2, BootHelper boothelper, Parameter parameter) {
    if (cipher1.logq - parameter.logp < parameter.logq) {
        boothelper.bootstrapping(cipher1, parameter.logq, parameter.logQ, parameter.logT);    
    }
    scheme.multAndEqual(cipher1, cipher2);
}

void fcnEncSqrtWithBoot(Ciphertext& outCipher, const Ciphertext& inCipher, Parameter parameter, long d, Scheme& scheme, BootHelper& boothelper) {
    long logp = parameter.logp;
    
    // cout << "=== start Sqrt ===" << endl;
    // cout << "initial logq = " << inCipher.logq << endl;
    // cout << "----------------" << endl;

    Ciphertext a = inCipher;
    Ciphertext b = inCipher;
    scheme.addConstAndEqual(b, -1.0, logp);

    // one iteration : 2logp 
    // total = 2 * d * logp
    Ciphertext dummy;
    for(int i = 0; i < d; i++) {

        cout << i << "/" << d - 1 << "th iteration" << endl;


        // ======================
        // Check the remaining logq first
        // ======================

        // cout << "Compare " << a.logq - 2 * logp << " vs " << parameter.logq << endl;
        if (a.logq - 2 * logp < parameter.logq) {
            cout << "bootstrapping in sqrt" << endl;
            boothelper.bootstrapping(a, parameter.logq, parameter.logQ, parameter.logT);
            boothelper.bootstrapping(b, parameter.logq, parameter.logQ, parameter.logT);
            // cout << "    after bootstrapping : logq = " << a.logq << ", " << b.logq << endl;
        }

        // make dummy = 1 - b / 2
        dummy = scheme.divByPo2(b, 1); // b - 1
        scheme.negateAndEqual(dummy); 
        // scheme.reScaleByAndEqual(dummy, logp); 
        scheme.addConstAndEqual(dummy, 1.0, logp);
        // dummy - 1

        // Update a
        // a <- a * (1 - b / 2)
        scheme.modDownToAndEqual(a, dummy.logq); // a - 1
        scheme.multAndEqual(a, dummy); 
        scheme.reScaleByAndEqual(a, logp); // a - logp + 1

        // make dummy = (b - 3) / 4
        dummy = scheme.addConst(b, -3.0, logp);
        // scheme.multByConstAndEqual(dummy, 0.25, logp);
        scheme.divByPo2AndEqual(dummy, 2); // dummy - 3
        // scheme.reScaleByAndEqual(dummy, logp); // dummy - logp

        //update b
        // b<- b * b * (b - 3) / 4
        // scheme.modDownByAndEqual(b, logp); // b - logp
        scheme.squareAndEqual(b);
        scheme.reScaleByAndEqual(b, logp); // b - logp
        scheme.modDownToAndEqual(dummy, b.logq);
        scheme.multAndEqual(b, dummy);
        scheme.reScaleByAndEqual(b, logp); // b - 2logp
        scheme.modDownToAndEqual(a, b.logq);

        // cout << "cipher.logq = " << a.logq << endl;
    }

    outCipher = a; 
    // cout << "final logq = " << outCipher.logq << endl;
    // cout << "=== End sqrt with boot ===" << endl;
}

void fcnEncMaxMinWithBoot(Ciphertext& maxCipher, Ciphertext& minCipher, Ciphertext& input1, Ciphertext& input2, Parameter parameter, long iter, Scheme& scheme, BootHelper& boothelper) {
    
    // cout << "=== start MaxMin ===" << endl;
    // cout << "initial logq = " << input1.logq << ", " << input2.logq << endl;
    // cout << "----------------" << endl;

    
    Ciphertext x = scheme.add(input1, input2);
    Ciphertext y = scheme.sub(input1, input2);
    scheme.divByPo2AndEqual(x, 1); // x - 1
    scheme.divByPo2AndEqual(y, 1); // y - 1

    // cout << "Compare " << x.logq - parameter.logp << " vs " << parameter.logq << endl;
    if (x.logq - parameter.logp < parameter.logq) {
        cout << " ---- Bootstrapping in MaxMin ---- " << endl;
        cout << x.logq - parameter.logp << " vs " << parameter.logq << endl;
        boothelper.bootstrapping(x, parameter.logq, parameter.logQ, parameter.logT);
        boothelper.bootstrapping(y, parameter.logq, parameter.logQ, parameter.logT);
    }
    
    scheme.squareAndEqual(y);
    scheme.reScaleByAndEqual(y, parameter.logp); // y - logp + 1

    Ciphertext sqrtCipher;
    // sqrtCipher - (2 * iter + 1) * logp + 1
    fcnEncSqrtWithBoot(sqrtCipher, y, parameter, iter, scheme, boothelper);

    if (x.logq < sqrtCipher.logq) {
        boothelper.bootstrapping(x, parameter.logq, parameter.logQ, parameter.logT);
    }
    scheme.modDownToAndEqual(x, sqrtCipher.logq);

    maxCipher = scheme.add(x, sqrtCipher);
    minCipher = scheme.sub(x, sqrtCipher);

    // cout << "----------------" << endl;
    // cout << "final logq = " << maxCipher.logq << ", " << minCipher.logq << endl;
    // cout << "=== End MaxMin ===" << endl;    
}

void fcnEncCompAndSwapWithBoot(Ciphertext& outCipher, const Ciphertext& inCipher, double* mask, long dist, Parameter parameter, long iter, Scheme& scheme, BootHelper& boothelper) {
    long logp = parameter.logp;
    Ciphertext a = inCipher;
    long n = a.n;

    // cout << "=== start CompAndSwap ===" << endl;
    // cout << "initial logq = " << inCipher.logq << endl;
    // cout << "----------------" << endl;

    // cout << "Compare " << a.logq - parameter.logp << " vs " << parameter.logq << endl;
    if (a.logq - parameter.logp < parameter.logq) {
        cout << " ---- Bootstrapping in CompAndSwap ---- " << endl;
        boothelper.bootstrapping(a, parameter.logq, parameter.logQ, parameter.logT);
    }
    
    Ciphertext maskCipher = scheme.encrypt(mask, n, logp, a.logq);
    Ciphertext dummy = scheme.mult(a, maskCipher);
    scheme.reScaleByAndEqual(dummy, logp); // dummy - logp
    scheme.modDownByAndEqual(a, logp);
    scheme.subAndEqual(a, dummy);
    scheme.rightRotateFastAndEqual(dummy, dist);

    Ciphertext max, min;
    fcnEncMaxMinWithBoot(max, min, a, dummy, parameter, iter, scheme, boothelper);
    
    scheme.leftRotateFastAndEqual(min, dist);
    scheme.addAndEqual(max, min);
    outCipher = max;

    // cout << "----------------" << endl;
    // cout << "final logq = " << outCipher.logq << endl;
    // cout << "=== End CompAndSwap ===" << endl;
}

long encBatcherOddEvenSortRec(Ciphertext& sortedCipher, const Ciphertext& inCipher, long logNum, long logJump, long loc, Parameter parameter, long iter, double** mask, Scheme& scheme, BootHelper& boothelper) {
    Ciphertext a = inCipher;
    Ciphertext dummy;
    if (logNum == 1) {
        cout << " <<<<<<< start "<< loc << "th CompAndSwap" << " <<<<<< " << endl;
        fcnEncCompAndSwapWithBoot(dummy, a, mask[loc], 1 << logJump, parameter, iter, scheme, boothelper);
        cout << " >>>>>>> end "<< loc << "th CompAndSwap" << " >>>>>> " << endl;
        a = dummy;
    } else {
        if (logJump == 0) {
            loc = encBatcherOddEvenSortRec(dummy, a, logNum - 1, logJump, loc, parameter, iter, mask, scheme, boothelper);
            a = dummy;
        }
        loc = encBatcherOddEvenSortRec(dummy, a, logNum - 1, logJump + 1, loc, parameter, iter, mask, scheme, boothelper);
        a = dummy;
        cout << " <<<<<<< start "<< loc << "th CompAndSwap" << " <<<<<< " << endl;
        fcnEncCompAndSwapWithBoot(dummy, a, mask[loc], 1 << logJump, parameter, iter, scheme, boothelper);
        cout << " >>>>>>> end "<< loc << "th CompAndSwap" << " >>>>>> " << endl;
        a = dummy;
    }
    sortedCipher = a;
    return loc + 1;
}

void encBatcherOddEvenSort(Ciphertext& sortedCipher, const Ciphertext& inCipher, Parameter parameter, long iter, Scheme& scheme, BootHelper& boothelper) {
    long log2n = parameter.log2n;
    double** mask;
    genAllMasking(log2n, mask);

    encBatcherOddEvenSortRec(sortedCipher, inCipher, log2n, 0, 0, parameter, iter, mask, scheme, boothelper);
}

void encBatcherOddEvenSortWithDecrypt(Ciphertext& sortedCipher, const Ciphertext& inCipher, Parameter parameter, long iter, Scheme& scheme, BootHelper& boothelper, SecretKey& secretKey) {
    long log2n = parameter.log2n;
    double** mask;
    genAllMasking(log2n, mask);

    encBatcherOddEvenSortRecWithDecrypt(sortedCipher, inCipher, log2n, 0, 0, parameter, iter, mask, scheme, boothelper, secretKey);
}


long encBatcherOddEvenSortRecWithDecrypt(Ciphertext& sortedCipher, const Ciphertext& inCipher, long logNum, long logJump, long loc, Parameter parameter, long iter, double** mask, Scheme& scheme, BootHelper& boothelper, SecretKey& secretKey) {
    Ciphertext a = inCipher;
    Ciphertext dummy;
    cout << "run encBatcherSort with loc = " << loc << endl;
   
    // cout << "call Rec(" << logNum << ", " << logJump << ")" << endl;
    if (logNum == 1) {
        // cout << "mask[" << loc << "] <- S(" << logNum << ", " << logJump << ")" << endl;
        // CompAndSwap()
        fcnEncCompAndSwapWithBootAndDecrypt(dummy, a, mask[loc], 1 << logJump, parameter, iter, scheme, boothelper, secretKey);

        // if (a.logq < 1500) {
        //     cout << "run Bootstrapping..." << endl;
        //     boothelper.bootstrapping(a, parameter.logq, parameter.logQ, parameter.logT);
        // }

        // fcnEncCompAndSwap(dummy, a, mask[loc], 1 << logJump, parameter.logp, iter, scheme);
        
        a = dummy;
        // compAndSwap(loc, 1 << logJump);
        // mask[loc] = genMaskingComp(1 << log2n, 1 << logJump);
    } else {
        if (logJump == 0) {
            loc = encBatcherOddEvenSortRecWithDecrypt(dummy, a, logNum - 1, logJump, loc, parameter, iter, mask, scheme, boothelper, secretKey);
            a = dummy;
            // loc = BatcherOddEvenSortRec(logNum - 1, logJump, loc);
        }
        loc = encBatcherOddEvenSortRecWithDecrypt(dummy, a, logNum - 1, logJump + 1, loc, parameter, iter, mask, scheme, boothelper, secretKey);
        a = dummy;
        // loc = BatcherOddEvenSortRec(logNum - 1, logJump + 1, loc);
        // cout << "mask[" << loc << "] <- S(" << logNum << ", " << logJump << ")" << endl;
        // mask[loc] = genMaskingMerge(1 << log2n, 1 << logNum, 1 << logJump);

        
        fcnEncCompAndSwapWithBootAndDecrypt(dummy, a, mask[loc], 1 << logJump, parameter, iter, scheme, boothelper, secretKey);

        // if (a.logq < 600) {
        //     cout << "run Bootstrapping..." << endl;
        //     boothelper.bootstrapping(a, parameter.logq, parameter.logQ, parameter.logT);
        // }
        // fcnEncCompAndSwap(dummy, a, mask[loc], 1 << logJump, parameter.logp, iter, scheme);
        
        a = dummy;
    }
    sortedCipher = a;
    return loc + 1;
}

void fcnEncCompAndSwapWithBootAndDecrypt(Ciphertext& outCipher, Ciphertext& inCipher, double* mask, long dist, Parameter parameter, long iter, Scheme& scheme, BootHelper& boothelper, SecretKey& secretKey) {
    long logp = parameter.logp;
    Ciphertext a = inCipher;
    long n = a.n;

    cout << "=== start CompAndSwap ===" << endl;
    cout << "inCipher.logq = " << inCipher.logq << endl;
    fcnDecryptAndPrint("inCipher", inCipher, scheme, secretKey);
    cout << "----------------" << endl;

    cout << "Compare " << a.logq - parameter.logp << " vs " << parameter.logq << endl;
    if (a.logq - 2 * parameter.logp < parameter.logq) {
        cout << " ---- Bootstrapping in CompAndSwap ---- " << endl;
        
        // scheme.multByConstAndEqual(a, 10000, logp);
        // scheme.reScaleByAndEqual(a, logp);
        
        boothelper.bootstrapping(a, parameter.logq, parameter.logQ, parameter.logT);
        
        // scheme.multByConstAndEqual(a, 0.0001, logp);
        // scheme.reScaleByAndEqual(a, logp);
    }
    
    Ciphertext maskCipher = scheme.encrypt(mask, n, logp, a.logq);
    Ciphertext dummy = scheme.mult(a, maskCipher);
    scheme.reScaleByAndEqual(dummy, logp); // dummy - logp
    scheme.modDownByAndEqual(a, logp);
    scheme.subAndEqual(a, dummy);
    scheme.rightRotateFastAndEqual(dummy, dist);
    
    Ciphertext max, min;
    fcnEncMaxMinWithBootAndDecrypt(max, min, a, dummy, parameter, iter, scheme, boothelper, secretKey);

    scheme.leftRotateFastAndEqual(min, dist);    
    scheme.addAndEqual(max, min);
    outCipher = max;
    
    cout << "outCipher.logq = " << outCipher.logq << endl;
    fcnDecryptAndPrint("outCipher", outCipher, scheme, secretKey);
    

    cout << "----------------" << endl;
    cout << "final logq = " << outCipher.logq << endl;
    cout << "=== End CompAndSwap ===" << endl;
}

void fcnEncMaxMinWithBootAndDecrypt(Ciphertext& maxCipher, Ciphertext& minCipher, Ciphertext& input1, Ciphertext& input2, Parameter parameter, long iter, Scheme& scheme, BootHelper& boothelper, SecretKey& secretKey) {
    
    cout << "=== start MaxMin ===" << endl;
    cout << "initial logq = " << input1.logq << ", " << input2.logq << endl;
    cout << "----------------" << endl;

    fcnDecryptAndPrintTwo("inputs", input1, input2, scheme, secretKey);
    
    Ciphertext x = scheme.add(input1, input2);
    Ciphertext y = scheme.sub(input1, input2);
    scheme.divByPo2AndEqual(x, 1); // x - 1
    scheme.divByPo2AndEqual(y, 1); // y - 1

    cout << "Compare " << x.logq - parameter.logp << " vs " << parameter.logq << endl;
    if (x.logq - 2 * parameter.logp < parameter.logq) {
        cout << " ---- Bootstrapping in MaxMin ---- " << endl;

        // scheme.multByConstAndEqual(x, 10000, parameter.logp);
        // scheme.reScaleByAndEqual(x, parameter.logp);
        // scheme.multByConstAndEqual(y, 10000, parameter.logp);
        // scheme.reScaleByAndEqual(y, parameter.logp);

        boothelper.bootstrapping(x, parameter.logq, parameter.logQ, parameter.logT);
        boothelper.bootstrapping(y, parameter.logq, parameter.logQ, parameter.logT);

        // scheme.multByConstAndEqual(x, 0.0001, parameter.logp);
        // scheme.reScaleByAndEqual(x, parameter.logp);
        // scheme.multByConstAndEqual(y, 0.0001, parameter.logp);
        // scheme.reScaleByAndEqual(y, parameter.logp);
    }
    
    scheme.squareAndEqual(y);
    scheme.reScaleByAndEqual(y, parameter.logp); // y - logp + 1

    Ciphertext sqrtCipher;
    // sqrtCipher - (2 * iter + 1) * logp + 1
    fcnEncSqrtWithBootWithDecrypt(sqrtCipher, y, parameter, iter, scheme, boothelper, secretKey);

    if (x.logq < sqrtCipher.logq) {
        cout << "x.logq = " << x.logq << " < sqrtCipher.logq = " << sqrtCipher.logq << endl;
        cout << "Bootstrapp x" << endl;
        boothelper.bootstrapping(x, parameter.logq, parameter.logQ, parameter.logT);
    }
    
    scheme.modDownToAndEqual(x, sqrtCipher.logq);

    maxCipher = scheme.add(x, sqrtCipher);
    minCipher = scheme.sub(x, sqrtCipher);

    
    cout << "outCipher.logq = " << maxCipher.logq << endl;
    fcnDecryptAndPrintTwo("MaxMin", maxCipher, minCipher, scheme, secretKey);

    cout << "----------------" << endl;
    cout << "final logq = " << maxCipher.logq << ", " << minCipher.logq << endl;
    cout << "=== End MaxMin ===" << endl;    
}


void fcnEncSqrtWithBootWithDecrypt(Ciphertext& outCipher, Ciphertext& inCipher, Parameter parameter, long d, Scheme& scheme, BootHelper& boothelper, SecretKey& secretKey) {
    long logp = parameter.logp;
    
    cout << "=== start Sqrt ===" << endl;
    cout << "initial logq = " << inCipher.logq << endl;
    fcnDecryptAndPrint("initial : ", inCipher, scheme, secretKey);
    cout << "----------------" << endl;

    Ciphertext a = inCipher;
    Ciphertext b = inCipher;
    scheme.addConstAndEqual(b, -1.0, logp);

    // one iteration : 2logp 
    // total = 2 * d * logp
    Ciphertext dummy;
    for(int i = 0; i < d; i++) {

        cout << i << "/" << d - 1 << "th iteration" << endl;


        // ======================
        // Check the remaining logq first
        // ======================

        cout << "Compare " << a.logq - 2 * logp << " vs " << parameter.logq << endl;
        if (a.logq - 2 * logp < parameter.logq) {

            // scheme.multByConstAndEqual(a, 10000, parameter.logp);
            // scheme.reScaleByAndEqual(a, parameter.logp);
            // scheme.multByConstAndEqual(b, 10000, parameter.logp);
            // scheme.reScaleByAndEqual(b, parameter.logp);

            boothelper.bootstrapping(a, parameter.logq, parameter.logQ, parameter.logT);
            boothelper.bootstrapping(b, parameter.logq, parameter.logQ, parameter.logT);

            // scheme.multByConstAndEqual(a, 0.0001, parameter.logp);
            // scheme.reScaleByAndEqual(a, parameter.logp);
            // scheme.multByConstAndEqual(b, 0.0001, parameter.logp);
            // scheme.reScaleByAndEqual(b, parameter.logp);
            
            cout << "    after bootstrapping : logq = " << a.logq << ", " << b.logq << endl;
        }

        // make dummy = 1 - b / 2
        dummy = scheme.divByPo2(b, 1); // b - 1
        scheme.negateAndEqual(dummy); 
        // scheme.reScaleByAndEqual(dummy, logp); 
        scheme.addConstAndEqual(dummy, 1.0, logp);
        // dummy - 1

        // Update a
        // a <- a * (1 - b / 2)
        scheme.modDownToAndEqual(a, dummy.logq); // a - 1
        scheme.multAndEqual(a, dummy); 
        scheme.reScaleByAndEqual(a, logp); // a - logp + 1

        // make dummy = (b - 3) / 4
        dummy = scheme.addConst(b, -3.0, logp);
        // scheme.multByConstAndEqual(dummy, 0.25, logp);
        scheme.divByPo2AndEqual(dummy, 2); // dummy - 3
        // scheme.reScaleByAndEqual(dummy, logp); // dummy - logp

        //update b
        // b<- b * b * (b - 3) / 4
        // scheme.modDownByAndEqual(b, logp); // b - logp
        scheme.squareAndEqual(b);
        scheme.reScaleByAndEqual(b, logp); // b - logp
        scheme.modDownToAndEqual(dummy, b.logq);
        scheme.multAndEqual(b, dummy);
        scheme.reScaleByAndEqual(b, logp); // b - 2logp
        scheme.modDownToAndEqual(a, b.logq);

        cout << "cipher.logq = " << a.logq << endl;
        fcnDecryptAndPrint("iter " + to_string(i) + " a :", a, scheme, secretKey);
    }

    outCipher = a; 
    cout << "final logq = " << outCipher.logq << endl;
    fcnDecryptAndPrint("output : ", outCipher, scheme, secretKey);
    cout << "=== End sqrt with boot ===" << endl;
}