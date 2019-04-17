#include "BootAlgorithm.h"

BootAlgo::BootAlgo(Parameter _param, long _iter, bool _increase) {
    param = _param;
    iter = _iter;
    increase = _increase;
}
void BootAlgo::approxSqrt(Ciphertext& cipher, BootScheme& scheme, BootHelper& bootHelper) {
    PrintUtils::nprint("start BootAlgo::sqrt", WANT_TO_PRINT);

    Ciphertext b = cipher;
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

        PrintUtils::nprint("After iter : logq = " + to_string(cipher.logq), WANT_TO_PRINT);
    }
    PrintUtils::nprint("end Sqrt, logq = " + to_string(cipher.logq), WANT_TO_PRINT);
}

void BootAlgo::minMax(Ciphertext& minCipher, Ciphertext& maxCipher, BootScheme& scheme, BootHelper& bootHelper) {
    PrintUtils::nprint("start minMax with logq = " + to_string(minCipher.logq) + ", " + to_string(maxCipher.logq), WANT_TO_PRINT);
    Ciphertext x = scheme.add(minCipher, maxCipher);
    Ciphertext y = scheme.sub(minCipher, maxCipher);
    scheme.divByPo2AndEqual(x, 1); // x - logp + 1
    scheme.divByPo2AndEqual(y, 1); // y - logp + 1
    
    scheme.squareAndEuqalWithBoot(y, bootHelper, param);
    scheme.reScaleByAndEqual(y, param.logp); // y - logp + 1

    // sqrtCipher - (2 * iter + 1) * logp + 1
    approxSqrt(y, scheme, bootHelper);

    // scheme.modDownToAndEqual(x, sqrtCipher.logq);
    scheme.modDownToAndEqualModified(x, y, bootHelper, param);

    maxCipher = scheme.add(x, y);
    minCipher = scheme.sub(x, y);
    PrintUtils::nprint("end minMax", WANT_TO_PRINT);
}

void BootAlgo::compAndSwap(Ciphertext& cipher, double* mask, long dist, BootScheme& scheme, Ring& ring, BootHelper& bootHelper) {
    PrintUtils::nprint("start compAndSwap with logq = " + to_string(cipher.logq), WANT_TO_PRINT);
    long n = cipher.n;
    ZZ* maskPoly = new ZZ[1 << param.logN];
    ring.encode(maskPoly, mask, n, param.logp);
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
        cout << "run compAndSwap with i = " << i << ", inc = " << increase << endl;
        compAndSwap(cipher, mask[i], 1 << (param.log2n - 1 - i), scheme, ring, bootHelper);
    }
    // compAndSwap(cipher, mask[0], 1 << (param.log2n - 1), scheme, ring, bootHelper);
}