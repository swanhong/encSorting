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

    double* mvec = EvaluatorUtils::randomRealArray(n, 0.25);
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
    for(int i = 0; i < 1 << param.log2n; i++) {
        cout << i << " : " << mvec[i] << " .. " << dvec[i].real() << ".." << mvec[i] - dvec[i].real() << endl;
    }
    PrintUtils::averageDifference(mvec, dvec, n);
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

    PlainSort plainSort;
    plainSort.runPlainSorting(mvec1, param.log2n, false);
    plainSort.runPlainSorting(mvec2, param.log2n, true);
	Ciphertext cipher1 = scheme.encrypt(mvec1, n, param.logp, param.logQ);
    Ciphertext cipher2 = scheme.encrypt(mvec2, n, param.logp, param.logQ);

    EncSorting encSorting(param, iter);
    
    // timeutils.start("EncSort");    
    // encSorting.runEncSorting(cipher1, scheme, ring, bootHelper, false);
    // encSorting.runEncSorting(cipher2, scheme, ring, bootHelper, true);
    // timeutils.stop("EncSort");

    Ciphertext* cipher = new Ciphertext[2];
    cipher[0] = cipher1;
    cipher[1] = cipher2;
    timeutils.start("Bitonic Merge");
    encSorting.bitonicMerge(cipher, 1, scheme, ring, bootHelper);
    timeutils.stop("Bitonic Merge"); 

    // // run PlainSort
    // CyclicArray ca(mvec, 1 << param.log2n);
    // PlainSort plainSort;
    // plainSort.runPlainSorting(ca, param.log2n, increase);
    // mvec = ca.getArray();

    // Print Result and Difference //	
	complex<double>* dvec1 = scheme.decrypt(secretKey, cipher[0]);
    complex<double>* dvec2 = scheme.decrypt(secretKey, cipher[1]);
    // PrintUtils::averageDifference(mvec, dvec, n);
    for(int i = 0; i < 1 << param.log2n; i++) {
        cout << i << " : " << mvec1[i] << ", " << dvec1[i].real() << endl;
    }
    cout << " === " << endl;
    for(int i = 0; i < 1 << param.log2n; i++) {
        cout << i << " : " << mvec2[i] << ", " << dvec2[i].real() << endl;
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
    CyclicArray* ca = new CyclicArray[n];
    for(int i = 0; i < num; i++) {
        // mvec[i] = EvaluatorUtils::randomRealArray(n);
        mvec[i] = new double[n];
        for(int j = 0; j < n; j++)
        {
            if (i % 2 == 0) {
                mvec[i][j] = (double) i * 0.25 + (double)j / 100 + 0.3;
            } else {
                mvec[i][j] =  (double) i * 0.25 + 0.51 - (double) j / 100 + 0.3;
            }
        }
        
        ca[i] = CyclicArray(mvec[i], 1 << param.log2n);
        plainSort.runPlainSorting(ca[i], param.log2n, (i % 2 == 0));
        ca[i].printAsVector();
        cipher[i] = scheme.encrypt(ca[i].getArray(), n, param.logp, param.logQ);
    }

    EncSorting encSorting(param, iter);
    timeutils.start("Bitonic Merge");
    encSorting.bitonicMerge(cipher, logNum, scheme, ring, bootHelper);
    timeutils.stop("Bitonic Merge"); 

    // run PlainMerge
    plainSort.bitonicMerge(ca, param.log2n, logNum);

    // Print Result and Difference //	
    for(int i = 0; i < num; i++) {
        complex<double>* dvec = scheme.decrypt(secretKey, cipher[i]);
        for(int j = 0; j < (1 << param.log2n); j++) {
            cout << j << " : " << ca[i].get(j) << ", " << dvec[j].real() << endl;
        }
        cout << "===" << endl;
    }
    
    
    // PrintUtils::averageDifference(mvec, dvec, n);
    
}