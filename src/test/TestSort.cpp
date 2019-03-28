#include "TestSort.h"

void TestSort::sort(Parameter param, long iter, bool increase) {
    srand(time(NULL));
	SetNumThreads(8);
    TimeUtils timeutils;
    
	PrintUtils::parameter(param, "TestSort");
    long n = 1 << param.log2n;
    long logp = param.logp;

    timeutils.start("TestSort KeyGen");
    Ring ring(param.logN, param.logQ);
    SecretKey secretKey(ring);
    BootScheme scheme(secretKey, ring);
    scheme.addConjKey(secretKey);
    scheme.addLeftRotKeys(secretKey);
    scheme.addRightRotKeys(secretKey);
    timeutils.stop("TestSort KeyGen");

	timeutils.start("Bootstrapping Helper construct");
	BootHelper bootHelper(param.log2n, param.radix, param.logc, scheme, ring, secretKey);
	timeutils.stop("Bootstrapping Helper construct");

    double* mvec = EvaluatorUtils::randomRealArray(n);
	Ciphertext cipher = scheme.encrypt(mvec, n, param.logp, param.logQ);

    timeutils.start("EncSort");
    EncSorting encSorting(param, iter);
    encSorting.runEncSorting(cipher, scheme, ring, bootHelper, increase);
    timeutils.stop("EncSort"); 

    // run PlainSort
    CyclicArray ca(mvec, 1 << param.log2n);
    PlainSort plainSort;
    plainSort.runPlainSorting(ca, param.log2n, increase);
    mvec = ca.getArray();

    // Print Result and Difference //	
	complex<double>* dvec = scheme.decrypt(secretKey, cipher);
    PrintUtils::averageDifference(mvec, dvec, n);
    // for(int i = 0; i < 1 << param.log2n; i++) {
    //     cout << i << " : " << mvec[i] << " .. " << dvec[i].real() << endl;
    // }
}

void TestSort::bitonicMerge(Parameter param, long iter) {
    srand(time(NULL));
	SetNumThreads(8);
    TimeUtils timeutils;
    
	PrintUtils::parameter(param, "Bitonic Merge");
    long n = 1 << param.log2n;
    long logp = param.logp;

    timeutils.start("TestSort KeyGen");
    Ring ring(param.logN, param.logQ);
    SecretKey secretKey(ring);
    BootScheme scheme(secretKey, ring);
    scheme.addConjKey(secretKey);
    scheme.addLeftRotKeys(secretKey);
    scheme.addRightRotKeys(secretKey);
    timeutils.stop("TestSort KeyGen");

	timeutils.start("Bootstrapping Helper construct");
	BootHelper bootHelper(param.log2n, param.radix, param.logc, scheme, ring, secretKey);
	timeutils.stop("Bootstrapping Helper construct");

    double* mvec1 = EvaluatorUtils::randomRealArray(n);
    double* mvec2 = EvaluatorUtils::randomRealArray(n);
	Ciphertext cipher1 = scheme.encrypt(mvec1, n, param.logp, param.logQ);
    Ciphertext cipher2 = scheme.encrypt(mvec2, n, param.logp, param.logQ);

    timeutils.start("EncSort");
    EncSorting encSorting(param, iter);
    encSorting.runEncSorting(cipher1, scheme, ring, bootHelper, false);
    encSorting.runEncSorting(cipher2, scheme, ring, bootHelper, true);
    timeutils.stop("EncSort");

    timeutils.start("Bitonic Merge");
    // encSorting.BitonicMerge(cipher1, cipher2, scheme, ring, bootHelper);
    timeutils.stop("Bitonic Merge"); 

    // // run PlainSort
    // CyclicArray ca(mvec, 1 << param.log2n);
    // PlainSort plainSort;
    // plainSort.runPlainSorting(ca, param.log2n, increase);
    // mvec = ca.getArray();

    // Print Result and Difference //	
	complex<double>* dvec1 = scheme.decrypt(secretKey, cipher1);
    complex<double>* dvec2 = scheme.decrypt(secretKey, cipher2);
    // PrintUtils::averageDifference(mvec, dvec, n);
    for(int i = 0; i < 1 << param.log2n; i++) {
        cout << i << " : " << dvec1[i].real() << endl;
    }
    cout << " === " << endl;
    for(int i = 0; i < 1 << param.log2n; i++) {
        cout << i << " : " << dvec2[i].real() << endl;
    }
}

void TestSort:: testMerge(Parameter param, long iter, long logNum) {
    srand(time(NULL));
	SetNumThreads(8);
    TimeUtils timeutils;
    long num = 1 << logNum;
    
	PrintUtils::parameter(param, "Merge");
    long n = 1 << param.log2n;
    long logp = param.logp;

    timeutils.start("KeyGen");
    Ring ring(param.logN, param.logQ);
    SecretKey secretKey(ring);
    BootScheme scheme(secretKey, ring);
    scheme.addConjKey(secretKey);
    scheme.addLeftRotKeys(secretKey);
    scheme.addRightRotKeys(secretKey);
    timeutils.stop("KeyGen");

	timeutils.start("Bootstrapping Helper construct");
	BootHelper bootHelper(param.log2n, param.radix, param.logc, scheme, ring, secretKey);
	timeutils.stop("Bootstrapping Helper construct");

    double ** mvec = new double*[num];
    Ciphertext* cipher = new Ciphertext[num];
    PlainSort plainSort;
    for(int i = 0; i < num; i++) {
        mvec[i] = EvaluatorUtils::randomRealArray(n);
        CyclicArray ca(mvec[i], n);
        plainSort.runPlainSorting(ca, param.log2n, (i % 2 == 0));
        ca.printAsVector();
        cipher[i] = scheme.encrypt(ca.getArray(), n, param.logp, param.logQ);
    }

    EncSorting encSorting(param, iter);
    timeutils.start("Bitonic Merge");
    encSorting.bitonicMerge(cipher, logNum, scheme, ring, bootHelper);
    timeutils.stop("Bitonic Merge"); 

    // // run PlainSort
    // CyclicArray ca(mvec, 1 << param.log2n);
    // PlainSort plainSort;
    // plainSort.runPlainSorting(ca, param.log2n, increase);
    // mvec = ca.getArray();

    // Print Result and Difference //	
    for(int i = 0; i < num; i++) {
        complex<double>* dvec = scheme.decrypt(secretKey, cipher[i]);
        for(int i = 0; i < 1 << param.log2n; i++) {
            cout << i << " : " << dvec[i].real() << endl;
        }
        cout << "===" << endl;
    }
    
    
    // PrintUtils::averageDifference(mvec, dvec, n);
    
}