

#include "TestAlgorithm.h"





// #include "../ImprovedHEAAN-master/HEAAN/app/common.h"



void TestAlgorithm::testMult(long logq, long logp, long logn) {

    srand(time(NULL));
    SetNumThreads(8);
    TimeUtils timeutils;
    Ring ring(logn, logq);
    SecretKey secretKey(ring);
    Scheme scheme(secretKey, ring);

    // **********************************

    long n = 1 << logn;
    complex<double>* mvec1 = EvaluatorUtils::randomComplexArray(n);
    complex<double>* mvec2 = EvaluatorUtils::randomComplexArray(n);
    complex<double>* mMult = new complex<double>[n];

    // **********************************
    
    for(long i = 0; i < n; i++) {
        mMult[i] = mvec1[i] * mvec2[i];
    }

    // **********************************
    
    Ciphertext cipher1 = scheme.encrypt(mvec1, n, logp, logq);
    Ciphertext cipher2 = scheme.encrypt(mvec2, n, logp, logq);

    // **********************************
    
    timeutils.start("Multiplication");
    scheme.multAndEqual(cipher1, cipher2);
    scheme.reScaleByAndEqual(cipher1, logp);
    timeutils.stop("Multiplication");

    // **********************************
    
    complex<double>* dmult = scheme.decrypt(secretKey, cipher1);
    StringUtils::compare(mMult, dmult, n, "mult");

    // **********************************
}

void TestAlgorithm::testPlainSort(long logn) {
    CyclicArray ca;
    ca.randomGen(1 << logn);

    cout << "Before sorting" << endl;
    ca.printAsVector();
    SortingAlgorithm sort(ca, logn);
    sort.BatcherOddEvenSort();
    cout << "After sorting" << endl;
    ca.printAsVector();
}

void TestAlgorithm::testMasking(long log2n) {
    double** mask;    
    cout << "gogo with " << log2n << endl;
    genAllMasking(log2n, mask);
    cout << "done" << endl;

    long length = 1 << log2n;
    
    for(int i = 0; i < (log2n + 1) * log2n / 2; i++){
        cout << "mask[" << i << "] = [";
        for(int j = 0; j < length; j++) {
            cout << mask[i][j] << ", ";
        }
        cout << "]" << endl;        
    }
}

void TestAlgorithm::testSqrt(Parameter parameter) {
	     // HE parameter //
    long logN = parameter.logN;
    long logQ = parameter.logQ;
    long logp = parameter.logp;
    long logc = parameter.logc;

    // Decomposition related parameter //
    long log2n = parameter.log2n;
    long radix = parameter.radix;

    // Bootstrapping parameter //
    long logq = parameter.logq;
    long logT = parameter.logT;

    long n = 1 << log2n;
	
	cout << "\n***************************" << endl;
    cout << "Test for encSqrt" << endl;
    cout << "logN = " << logN << ", logQ = " << logQ << ", logp = " << logp << ", logc = " << logc << endl;
    cout << "slots = " << n << ", radix = " << radix << ", logq = " << logq << ", logT = " << logT << endl;
    cout << "***************************" << endl;
    cout << endl;


    TimeUtils timeutils;
    timeutils.start("KeyGen");
    Ring ring(logN, logQ);
    SecretKey secretKey(ring);
    Scheme scheme(secretKey, ring);
    scheme.addConjKey(secretKey);
    scheme.addLeftRotKeys(secretKey);
    timeutils.stop("KeyGen");

    double* mvec = new double[n];
    for(int i = 0; i < n; i++) {
        mvec[i] = 0.00001 + (double) i / n;
    }
    
	Ciphertext cipher = scheme.encrypt(mvec, n, logp, logQ);

	Ciphertext cipher2;

    timeutils.start("Sqrt");
    fcnEncComputeSqrt(cipher2, cipher, logp, 5, scheme);
       
    timeutils.stop("Sqrt");

	complex<double>* dmult = scheme.decrypt(secretKey, cipher2);
    // StringUtils::compare(mvec, dmult, n, "mult");
    for(int i = 0; i < n; i++) {
        cout << i << " : " << mvec[i] << ", " << dmult[i].real() << endl;
    }
    

}

void TestAlgorithm::testMaxMin(Parameter parameter) {
         // HE parameter //
    long logN = parameter.logN;
    long logQ = parameter.logQ;
    long logp = parameter.logp;
    long logc = parameter.logc;

    // Decomposition related parameter //
    long log2n = parameter.log2n;
    long radix = parameter.radix;

    // Bootstrapping parameter //
    long logq = parameter.logq;
    long logT = parameter.logT;

    long n = 1 << log2n;
	
	cout << "\n***************************" << endl;
    cout << "Test for encMaxMin" << endl;
    cout << "logN = " << logN << ", logQ = " << logQ << ", logp = " << logp << ", logc = " << logc << endl;
    cout << "slots = " << n << ", radix = " << radix << ", logq = " << logq << ", logT = " << logT << endl;
    cout << "***************************" << endl;
    cout << endl;


    TimeUtils timeutils;
    timeutils.start("KeyGen");
    Ring ring(logN, logQ);
    SecretKey secretKey(ring);
    Scheme scheme(secretKey, ring);
    scheme.addConjKey(secretKey);
    scheme.addLeftRotKeys(secretKey);
    timeutils.stop("KeyGen");

    double* mvec1 = new double[n];
    for(int i = 0; i < n; i++) {
        mvec1[i] = 0.00001 + (double) i / n;
    }

    double* mvec2 = new double[n];
    for(int i = 0; i < n; i++) {
        mvec2[i] = 1.0 - (double) i / n;
    }
    
	Ciphertext cipher1 = scheme.encrypt(mvec1, n, logp, logQ);
    Ciphertext cipher2 = scheme.encrypt(mvec2, n, logp, logQ);
    
    Ciphertext minCipher, maxCipher;
    fcnEncMaxMin(maxCipher, minCipher, cipher1, cipher2, logp, 5, scheme);
    
    complex<double>* dvec1 = scheme.decrypt(secretKey, maxCipher);
    complex<double>* dvec2 = scheme.decrypt(secretKey, minCipher);
    // StringUtils::compare(mvec, dmult, n, "mult");
    for(int i = 0; i < n; i++) {
        cout << i << " : " << 
            mvec1[i] << ", " << mvec2[i] << " /// -> /// " <<
            dvec1[i].real() << ", " << dvec2[i].real() << endl;
    }
}

void TestAlgorithm::testEncCompAndSwap(Parameter parameter) {
         // HE parameter //
    long logN = parameter.logN;
    long logQ = parameter.logQ;
    long logp = parameter.logp;
    long logc = parameter.logc;

    // Decomposition related parameter //
    long log2n = parameter.log2n;
    long radix = parameter.radix;

    // Bootstrapping parameter //
    long logq = parameter.logq;
    long logT = parameter.logT;

    long n = 1 << log2n;
	
	cout << "\n***************************" << endl;
    cout << "Test for encCompAndSwap" << endl;
    cout << "logN = " << logN << ", logQ = " << logQ << ", logp = " << logp << ", logc = " << logc << endl;
    cout << "slots = " << n << ", radix = " << radix << ", logq = " << logq << ", logT = " << logT << endl;
    cout << "***************************" << endl;
    cout << endl;


    TimeUtils timeutils;
    timeutils.start("KeyGen");
    Ring ring(logN, logQ);
    SecretKey secretKey(ring);
    Scheme scheme(secretKey, ring);
    scheme.addConjKey(secretKey);
    scheme.addLeftRotKeys(secretKey);
    scheme.addRightRotKeys(secretKey);
    timeutils.stop("KeyGen");

    double* mvec = new double[n];
    for(int i = 0; i < n; i++) {
        mvec[n - 1 - i] = 0.00001 + (double) i / n;
    }
    
	Ciphertext cipher = scheme.encrypt(mvec, n, logp, logQ);

	Ciphertext cipher2;

    double* mask = genMaskingComp(n, 1);
    
    timeutils.start("CompAndSwap");
    fcnEncCompAndSwap(cipher2, cipher, mask, 1, logp, 7, scheme);
    timeutils.stop("CompAndSwap");

	complex<double>* dmult = scheme.decrypt(secretKey, cipher2);
    // StringUtils::compare(mvec, dmult, n, "mult");
    for(int i = 0; i < n; i++) {
        cout << i << " : " << mvec[i] << ", " << dmult[i].real() << endl;
    }
}

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

void fcnDecryptAndPrint(string str, Ciphertext cipher, Scheme scheme, SecretKey secretKey) {
    complex<double>* dvec = scheme.decrypt(secretKey, cipher);
    cout << "==== " << str << " ====" << endl;
    for(int i = 0; i < cipher.n; i++) {
        cout << dvec[i].real() << endl;
    }
}

