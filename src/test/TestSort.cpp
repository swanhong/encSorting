#include "TestSort.h"

void TestSort::sort(Parameter param, long iter) {
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
    Ciphertext sortedCipher;

    timeutils.start("EncSort");
    EncSorting encSorting;
    encSorting.runSorting(sortedCipher, cipher, param, iter, scheme, ring, bootHelper);
    timeutils.stop("EncSort"); 

    // run PlainSort
    CyclicArray ca(mvec, 1 << param.log2n);
    PlainSort plainSort;
    plainSort.runSorting(ca, param.log2n);
    mvec = ca.getArray();

    // Print Result and Difference //	
	complex<double>* dvec = scheme.decrypt(secretKey, sortedCipher);
    PrintUtils::averageDifference(mvec, dvec, n);
}