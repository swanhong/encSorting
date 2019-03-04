#include "../HEAAN/src/HEAAN.h"
#include "../bootsrc/new_bootstrapping.h"
#include "Parameter.h"

void multWithBootAndEqual(Scheme& scheme, Ciphertext& cipher1, Ciphertext& cipher2, BootHelper boothelper, Parameter parameter);

void fcnEncSqrtwithBoot(Ciphertext& outCipher, const Ciphertext& inCipher, Parameter parameter, long d, Scheme scheme, BootHelper boothelper);

void fcnEncMaxMinWithBoot(Ciphertext& maxCipher, Ciphertext& minCipher, Ciphertext& input1, Ciphertext& input2, Parameter parameter, long iter, Scheme scheme, BootHelper boothelper);

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

void fcnEncSqrtwithBoot(Ciphertext& outCipher, const Ciphertext& inCipher, Parameter parameter, long d, Scheme scheme, BootHelper boothelper) {
    long logp = parameter.logp;
    
    cout << "start Sqrt" << endl;
    cout << "initial logq = " << inCipher.logq << endl;
    cout << "----------------" << endl;

    Ciphertext a = inCipher;
    Ciphertext b = inCipher;
    scheme.addConstAndEqual(b, -1.0, logp);

    // one iteration : 2logp 
    // total = 2 * d * logp
    Ciphertext dummy;
    for(int i = 0; i < d; i++) {

        std::cout << i << "/" << d - 1 << "th iteration" << '\n';


        // ======================
        // Check the remaining logq first
        // ======================

        cout << "Compare " << a.logq - 2 * logp << " vs " << parameter.logq << endl;
        if (a.logq - 2 * logp < parameter.logq) {
            boothelper.bootstrapping(a, parameter.logq, parameter.logQ, parameter.logT);
            boothelper.bootstrapping(b, parameter.logq, parameter.logQ, parameter.logT);
        }

        
        // make dummy = 1 - b / 2
        dummy = scheme.multByConst(b, -0.5, logp);
        scheme.reScaleByAndEqual(dummy, logp); 
        scheme.addConstAndEqual(dummy, 1.0, logp);
        // dummy - logp

        // Update a
        // a <- a * (1 - b / 2)
        scheme.modDownByAndEqual(a, logp); // a - logp
        scheme.multAndEqual(a, dummy); 
        scheme.reScaleByAndEqual(a, logp); // a - 2logp

        // make dummy = (b - 3) / 4
        dummy = scheme.addConst(b, -3.0, logp);
        scheme.multByConstAndEqual(dummy, 0.25, logp);
        scheme.reScaleByAndEqual(dummy, logp); // dummy - logp

        //update b
        // b<- b * b * (b - 3) / 4
        // scheme.modDownByAndEqual(b, logp); // b - logp
        scheme.squareAndEqual(b);
        scheme.reScaleByAndEqual(b, logp); // b - logp
        scheme.multAndEqual(b, dummy);
        scheme.reScaleByAndEqual(b, logp); // b - 2logp

        cout << "a.logq = " << a.logq << endl;
        cout << "b.logq = " << b.logq << endl;
    }

    outCipher = a; 
    cout << "End sqrt with boot" << endl;
}

void fcnEncMaxMinWithBoot(Ciphertext& maxCipher, Ciphertext& minCipher, Ciphertext& input1, Ciphertext& input2, Parameter parameter, long iter, Scheme scheme, BootHelper boothelper) {
    Ciphertext x = scheme.add(input1, input2);
    Ciphertext y = scheme.sub(input1, input2);
    scheme.divByPo2AndEqual(x, 1); // x - 1
    scheme.divByPo2AndEqual(y, 1); // y - 1

    cout << "x.logq = " << x.logq << endl;
    
    scheme.squareAndEqual(y);
    scheme.reScaleByAndEqual(y, parameter.logp); // y - logp + 1

    Ciphertext sqrtCipher;
    // sqrtCipher - (2 * iter + 1) * logp + 1
    fcnEncSqrtwithBoot(sqrtCipher, y, parameter, iter, scheme, boothelper);

    scheme.modDownToAndEqual(x, sqrtCipher.logq);

    cout << "asdf " << endl;
    maxCipher = scheme.add(x, sqrtCipher);
    cout << "asdf " << endl;
    minCipher = scheme.sub(x, sqrtCipher);
    cout << "asdf " << endl;
}