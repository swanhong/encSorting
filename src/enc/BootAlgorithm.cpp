#include "BootAlgorithm.h"


void BootAlgo::approxSqrt(Ciphertext& outCipher, const Ciphertext& inCipher, Parameter param, long iter, BootScheme& scheme, BootHelper& bootHelper) {
    PrintUtils::nprint("start BootAlgo::sqrt", WANT_TO_PRINT);

    Ciphertext a = inCipher;
    Ciphertext b = inCipher;
    long logp = param.logp;    
    scheme.addConstAndEqual(b, -1.0, logp);

    // one iteration : 2logp 
    // total = 2 * d * logp
    Ciphertext dummy;
    for(int i = 0; i < iter; i++) {
        PrintUtils::nprint(to_string(i) + "/" + to_string(iter - 1) + "th iteration", WANT_TO_PRINT);


        // make dummy = 1 - b / 2
        dummy = scheme.divByPo2(b, 1); // b - 1
        scheme.negateAndEqual(dummy); 
        scheme.addConstAndEqual(dummy, 1.0, logp); // dummy - 1

        // Update a
        // a <- a * (1 - b / 2)
        scheme.modDownToAndEqualModified(a, dummy, bootHelper, param); // a - 1
        scheme.multAndEqualWithBoot(a, dummy, bootHelper, param);
        scheme.reScaleByAndEqual(a, logp); // a - logp + 1

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
        scheme.modDownToAndEqualModified(a, b, bootHelper, param);

        PrintUtils::nprint("After iter : logq = " + to_string(a.logq), WANT_TO_PRINT);
    }

    outCipher = a; 
    cout << "end Sqrt, logq = " << outCipher.logq << endl;
}

void BootAlgo::minMax(Ciphertext& minCipher, Ciphertext& maxCipher, Ciphertext& input1, Ciphertext& input2, long iter, Parameter& param, BootScheme& scheme, BootHelper& bootHelper) {
    cout << "start minMax with logq = " << input1.logq << ", " << input2.logq << endl;
    Ciphertext x = scheme.add(input1, input2);
    Ciphertext y = scheme.sub(input1, input2);
    scheme.divByPo2AndEqual(x, 1); // x - logp + 1
    scheme.divByPo2AndEqual(y, 1); // y - logp + 1
    
    scheme.squareAndEuqalWithBoot(y, bootHelper, param);
    scheme.reScaleByAndEqual(y, param.logp); // y - logp + 1

    Ciphertext sqrtCipher;
    // sqrtCipher - (2 * iter + 1) * logp + 1
    approxSqrt(sqrtCipher, y, param, iter, scheme, bootHelper);

    // scheme.modDownToAndEqual(x, sqrtCipher.logq);
    scheme.modDownToAndEqualModified(x, sqrtCipher, bootHelper, param);

    maxCipher = scheme.add(x, sqrtCipher);
    minCipher = scheme.sub(x, sqrtCipher);
    cout << "end minMax" << endl;
}

void BootAlgo::compAndSwap(Ciphertext& outCipher, const Ciphertext& inCipher, double* mask, long dist, long iter, Parameter& param, BootScheme& scheme, Ring& ring, BootHelper& bootHelper) {
    cout << "start compAndSwap with logq = " << inCipher.logq << endl;
    Ciphertext a = inCipher;
    long n = a.n;
    ZZ* maskPoly = new ZZ[1 << param.logN];
    ring.encode(maskPoly, mask, n, param.logp);
    // Ciphertext maskCipher = scheme.encrypt(mask, n, param.logp, param.logQ);
    Ciphertext dummy = a;
    scheme.multByPolyAndEqualWithBoot(dummy, maskPoly, bootHelper, param);
    scheme.reScaleByAndEqual(dummy, param.logp);
    scheme.modDownToAndEqualModified(a, dummy, bootHelper, param);
    scheme.subAndEqual(a, dummy);
    scheme.rightRotateFastAndEqual(dummy, dist);
    Ciphertext min, max;
    minMax(min, max, a, dummy, iter, param, scheme,  bootHelper);
    scheme.leftRotateFastAndEqual(min, dist);
    scheme.addAndEqual(max, min);

    outCipher = max;
    cout << "end compAndSwap" << endl;
}