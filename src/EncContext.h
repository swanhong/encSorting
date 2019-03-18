#include "../HEAAN/src/HEAAN.h"

void fcnDecryptAndPrint(string str, Ciphertext& cipher, Scheme& scheme, SecretKey& secretKey);
void fcnDecryptAndPrintTwo(string str, Ciphertext& cipher1, Ciphertext& cipher2, Scheme& scheme, SecretKey& secretKey);

void fcnEncComputeSqrt(Ciphertext& outCipher, const Ciphertext& inCipher, long logp, long d, Scheme scheme);

void fcnEncMaxMin(Ciphertext& maxCipher, Ciphertext& minCipher, Ciphertext& input1, Ciphertext& input2, long logp, long iter, Scheme scheme);

void fcnGenEncAllMasking(long log2n, double**& mask);

long fcnGenEncMaskingRec(long log2n, long logNum, long logJump, double**& mask, long loc);

double* fcnGenEncMaskingComp(long length, long jump);

double* fcnGenEncMaskingMerge(long length, long num, long jump);

void fcnEncCompAndSwap(Ciphertext& outCipher, const Ciphertext& inCipher, double* mask, long dist, long logp, long iter, Scheme scheme);


void fcnEncComputeSqrt(Ciphertext& outCipher, const Ciphertext& inCipher, long logp, long d, Scheme scheme) {
    cout << "start Sqrt" << endl;
    Ciphertext a = inCipher;
    Ciphertext b = inCipher;
    scheme.addConstAndEqual(b, -1.0, logp);

    // one iteration : 2logp 
    // total = 2 * d * logp
    Ciphertext dummy;
    for(int i = 0; i < d; i++) {
        std::cout << i << "/" << d - 1 << "th iteration" << '\n';
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

        cout << "a.logq = " << a.logq << endl;
        cout << "b.logq = " << b.logq << endl;
    }

    outCipher = a;    
}
 
void fcnEncMaxMin(Ciphertext& maxCipher, Ciphertext& minCipher, Ciphertext& input1, Ciphertext& input2, long logp, long iter, Scheme scheme) {
    Ciphertext x = scheme.add(input1, input2);
    Ciphertext y = scheme.sub(input1, input2);
    scheme.divByPo2AndEqual(x, 1); // x - logp + 1
    scheme.divByPo2AndEqual(y, 1); // y - logp + 1

    cout << "x.logq = " << x.logq << endl;
    
    scheme.squareAndEqual(y);
    scheme.reScaleByAndEqual(y, logp); // y - logp + 1

    Ciphertext sqrtCipher;
    // sqrtCipher - (2 * iter + 1) * logp + 1
    fcnEncComputeSqrt(sqrtCipher, y, logp, iter, scheme);

    scheme.modDownToAndEqual(x, sqrtCipher.logq);

    maxCipher = scheme.add(x, sqrtCipher);
    minCipher = scheme.sub(x, sqrtCipher);
}

void fcnGenAllEncMasking(long log2n, double**& mask) {
    int maskNum = (log2n + 1) * log2n / 2;
    mask = new double*[maskNum];
    fcnGenEncMaskingRec(log2n, log2n, 0, mask, 0);
}

long fcnGenEncMaskingRec(long log2n, long logNum, long logJump, double**& mask, long loc) {

    // cout << "call Rec(" << log2n << ", " << logNum << ", " << logJump << ")" << endl;
    if (logNum == 1) {
        // cout << "mask[" << loc << "] <- S(" << logNum << ", " << logJump << ")" << endl;
        mask[loc] = fcnGenEncMaskingComp(1 << log2n, 1 << logJump);
        return loc + 1;
    } else {
        if (logJump == 0) {
            loc = fcnGenEncMaskingRec(log2n, logNum - 1, logJump, mask, loc);
        }
        loc = fcnGenEncMaskingRec(log2n, logNum - 1, logJump + 1, mask, loc);
        // cout << "mask[" << loc << "] <- S(" << logNum << ", " << logJump << ")" << endl;
        mask[loc] = fcnGenEncMaskingMerge(1 << log2n, 1 << logNum, 1 << logJump);
        return loc + 1;
    }
}

double* fcnGenEncMaskingComp(long length, long jump) {   
    double* mask = new double[length];
    for(int i = 0; i < length; i++) {
        mask[i] = 0;
    }
    
    long repeat = length / (jump * 2);
    for(int i = 0; i < repeat; i++) {
        for(int j = 0; j < jump; j++) {
            mask[i * jump * 2 + j] = 1;
        }
    }
    return mask;
}

double* fcnGenEncMaskingMerge(long length, long num, long jump) {
    double* mask = new double[length];
    for(int i = 0; i < length; i++) {
        mask[i] = 0;
    }
    long repeat = length / (jump * num);
    for(int i = 0; i < repeat; i++) {
        for(int j = 0; j < jump; j++) {
            for(int k = 0; k < num / 2 - 1; k++) {
                mask[i * jump * num + (2 * k + 1) * jump + j] = 1;    
            }
            
            
        }        
    }
    return mask;
}

void fcnEncCompAndSwap(Ciphertext& outCipher, const Ciphertext& inCipher, double* mask, long dist, long logp, long iter, Scheme scheme) {
    Ciphertext a = inCipher;
    long n = a.n;
    long logQ = a.logq;
    Ciphertext maskCipher = scheme.encrypt(mask, n, logp, logQ);
    Ciphertext dummy = scheme.mult(a, maskCipher);
    scheme.reScaleByAndEqual(dummy, logp);
    scheme.modDownByAndEqual(a, logp);
    scheme.subAndEqual(a, dummy);
    scheme.rightRotateFastAndEqual(dummy, dist);
    Ciphertext max, min;
    fcnEncMaxMin(max, min, a, dummy, logp, iter, scheme);
    scheme.leftRotateFastAndEqual(min, dist);
    scheme.addAndEqual(max, min);

    outCipher = max;
}

void fcnDecryptAndPrint(string str, Ciphertext& cipher, Scheme& scheme, SecretKey& secretKey) {
    complex<double>* dvec = scheme.decrypt(secretKey, cipher);
    cout << "==== " << str << " ====" << endl;

    long POWER_OF_TEN = 1;
    while(POWER_OF_TEN * 100 < cipher.n){
        POWER_OF_TEN *= 10;
    }
    
    cout << "num_slots = " << cipher.n << endl;
    cout << "POWER_OF_TEN = " << POWER_OF_TEN << endl;

    for(int i = 0; i < cipher.n; i++) {
        if (cipher.n > 100) {
            if (i % POWER_OF_TEN == 0) {
                cout << i << " : " << dvec[i] << endl;
            }        
        } else {
            cout << i << " : " << dvec[i] << endl;
        }      
    }
}

void fcnDecryptAndPrintTwo(string str, Ciphertext& cipher1, Ciphertext& cipher2, Scheme& scheme, SecretKey& secretKey) {
    complex<double>* dvec1 = scheme.decrypt(secretKey, cipher1);
    complex<double>* dvec2 = scheme.decrypt(secretKey, cipher2);
    
    cout << "==== " << str << " ====" << endl;
    for(int i = 0; i < cipher1.n; i++) {
        if (i % 100 == 0) {
            cout << i << " : " << dvec1[i].real() << ", " << dvec2[i].real() << endl;
        }        
    }
}