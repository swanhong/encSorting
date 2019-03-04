#include "EncContext.h"
#include "CyclicArray.h"
#include "SortingAlrogithm.h"
#include "Parameter.h"

class TestAlgorithm {
public:
    static void testBootstrapping(Parameter parameter);

    static void testOpreationsAndBoot(Parameter parameter);

    static void testMult(long logq, long logp, long logn);

    // static void testCompAndSwap(long length);

    static void testPlainSort(long logn);

    static void testMasking(long log2n);

    static void testSqrt(Parameter parameter);

    static void testMaxMin(Parameter parameter);

    static void testEncCompAndSwap(Parameter parameter);
};



// int fcnEncBatcherOddEvenSortRec(Ciphertext& cipher, long num, long logJump, long loc, double** mask, long iter, Scheme scheme, BootHelper boothelper, Parameter parameter);

// void fcnEncBatcherOddEvenSort(Ciphertext& sortedCipher, const Ciphertext& inCipher, double** mask, long iter, Scheme scheme, BootHelper boothelper, Parameter parameter);

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