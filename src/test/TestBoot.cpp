#include "TestBoot.h"

void TestBoot::bootstrapping(Parameter parameter) {
	PrintUtils::parameter(parameter, "Improved Bootstrapping");
	
	long n = 1 << parameter.log2n;   
	
	TimeUtils timeutils;
	timeutils.start("KeyGen");
	Ring ring(parameter.logN, parameter.logQ);
	SecretKey secretKey(ring);
	Scheme scheme(secretKey, ring);
	scheme.addConjKey(secretKey);
	scheme.addLeftRotKeys(secretKey);
	
	srand(time(NULL));
	SetNumThreads(8);
	timeutils.stop("KeyGen");

	timeutils.start("Bootstrapping Helper construct");
	BootHelper boothelper(parameter.log2n, parameter.radix, parameter.logc, scheme, ring, secretKey);
	timeutils.stop("Bootstrapping Helper construct");

    complex<double>* mvec = EvaluatorUtils::randomComplexArray(n);

    Ciphertext cipher1 = scheme.encrypt(mvec, n, parameter.logp, parameter.logQ);
    Ciphertext cipher2 = scheme.encrypt(mvec, n, parameter.logp, parameter.logQ);
    
    timeutils.start("Improved bootstrapping");
    boothelper.bootstrapping(cipher1, parameter.logq, parameter.logQ, 4, 4); 
    timeutils.stop("Improved bootstrapping");
    cout << "* consumed logQ = " << parameter.logQ - cipher1.logq << endl;
    complex<double>* dvec1 = scheme.decrypt(secretKey, cipher1);
    cout << "log2(avg of error) = " << diff(mvec, dvec1, n) << endl;
    
    timeutils.start("Improved cos bootstrappingg");
    boothelper.bootstrapping_cos(cipher2, parameter.logq, parameter.logQ, 8);
    timeutils.stop("Improved cos bootstrappingg");
    cout << "* consumed logQ = " << parameter.logQ - cipher2.logq << endl;
    complex<double>* dvec2 = scheme.decrypt(secretKey, cipher2);
    cout << "log2(avg of error) = " << diff(mvec, dvec2, n) << endl;
    
    // Print Result and Difference //
    
    // if(mvec != NULL) delete[] mvec;
    // if(dvec != NULL) delete[] dvec;
		
	return;
}

void TestBoot::approxSqrt(Parameter parameter, long iter) {
	long n = 1 << parameter.log2n;

	PrintUtils::parameter(parameter, "TestBoot::approxSqrt");
	
	TimeUtils timeutils;
	timeutils.start("KeyGen");
	Ring ring(parameter.logN, parameter.logQ);
	SecretKey secretKey(ring);
	BootScheme scheme(secretKey, ring);
	scheme.addConjKey(secretKey);
	scheme.addLeftRotKeys(secretKey);
	timeutils.stop("KeyGen");

	timeutils.start("Bootstrapping Helper construct");
	BootHelper bootHelper(parameter.log2n, parameter.radix, parameter.logc, scheme, ring, secretKey);
	timeutils.stop("Bootstrapping Helper construct");

	double* mvec = EvaluatorUtils::randomRealArray(n);
	Ciphertext cipher = scheme.encrypt(mvec, n, parameter.logp, parameter.logQ);

    
    timeutils.start("Sqrt");
    BootAlgo bootAlgo(parameter, iter);
	bootAlgo.approxSqrt(cipher, scheme, bootHelper);
    timeutils.stop("Sqrt");

    // Print Result and Difference //	
	complex<double>* dvec = scheme.decrypt(secretKey, cipher);
    for(int i = 0; i < n; i++) {
        mvec[i] = sqrt(mvec[i]);
    }
    PrintUtils::averageDifference(mvec, dvec, n);
}

void TestBoot::minMax(Parameter param, long iter) {
    PrintUtils::parameter(param, "TestBoot::minMax");
	long n = 1 << param.log2n;

    TimeUtils timeutils;
    timeutils.start("KeyGen");
    Ring ring(param.logN, param.logQ);
    SecretKey secretKey(ring);
    BootScheme scheme(secretKey, ring);
    scheme.addConjKey(secretKey);
    scheme.addLeftRotKeys(secretKey);
	srand(time(NULL));
	SetNumThreads(8);
    timeutils.stop("KeyGen");

	timeutils.start("Bootstrapping Helper construct");
	BootHelper boothelper(param.log2n, param.radix, param.logc, scheme, ring, secretKey);
	timeutils.stop("Bootstrapping Helper construct");

	double* mvec1 = EvaluatorUtils::randomRealArray(n);
    double* mvec2 = EvaluatorUtils::randomRealArray(n);
    
	Ciphertext cipher1 = scheme.encrypt(mvec1, n, param.logp, param.logQ);
    Ciphertext cipher2 = scheme.encrypt(mvec2, n, param.logp, param.logQ);
    
	timeutils.start("minMax");
    BootAlgo bootAlgo(param, iter);
    bootAlgo.minMax(cipher1, cipher2, scheme, boothelper);
    timeutils.stop("minMax");

	complex<double>* dmax = scheme.decrypt(secretKey, cipher2);
    complex<double>* dmin = scheme.decrypt(secretKey, cipher1);
    
    double* max = new double[n];
    for(int i = 0; i < n; i++) {
        if(mvec1[i] > mvec2[i]) {
            max[i] = mvec1[i];
        } else {
            max[i] = mvec2[i];
        }   
    }
    PrintUtils::averageDifference(max, dmax, n);
}

void TestBoot::compAndSwap(Parameter param, long iter) {
    srand(time(NULL));
	SetNumThreads(8);
    
    long n = 1 << param.log2n;
	
    PrintUtils::parameter(param, "TestBoot::compAndSwap");

    TimeUtils timeutils;
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


    double* mvec = EvaluatorUtils::randomRealArray(n);

	Ciphertext cipher = scheme.encrypt(mvec, n, param.logp, param.logQ);

    MaskingGenerator mg(param.log2n);
    double** mask = mg.getBitonicMergeMasking();

    BootAlgo bootAlgo(param, iter);
    PlainSort plainSort;
    for(int i = 0; i < param.log2n; i++)
    {
        timeutils.start("CompAndSwap");
        bootAlgo.compAndSwap(cipher, mask[i], 1 << (param.log2n - 1 - i), scheme, ring, bootHelper);	
        
        timeutils.stop("CompAndSwap");
        
        CyclicArray ca(mvec, n);
        plainSort.compAndSwap(ca, mask, i, 1 << (param.log2n - 1 - i), true);
        mvec = ca.getArray();
        // for(int i = 0; i < n / 2; i++) {
        //     if(mvec[i] > mvec[i + n / 2]) {
        //         double x = mvec[i];
        //         mvec[i] = mvec[i + n / 2];
        //         mvec[i + n / 2] = x;
        //     }
        // }

        complex<double>* dvec = scheme.decrypt(secretKey, cipher);
        PrintUtils::printArrays(mvec, dvec, n);
        PrintUtils::averageDifference(mvec, dvec, n);
    }
}

void TestBoot::compAndSwapWithBoot(Parameter param, long iter) {
    srand(time(NULL));
	SetNumThreads(8);
    
    long n = 1 << param.log2n;
	
    PrintUtils::parameter(param, "TestBoot::compAndSwap");

    TimeUtils timeutils;
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


    double* mvec = EvaluatorUtils::randomRealArray(n, 0.25);

	Ciphertext cipher = scheme.encrypt(mvec, n, param.logp, param.logQ);

    MaskingGenerator mg(param.log2n);
    double** mask = mg.getBitonicMergeMasking();

    BootAlgo bootAlgo(param, iter);
    PlainSort plainSort;
    for(int i = 0; i < param.log2n; i++)
    {
        timeutils.start("CompAndSwap");
        bootAlgo.compAndSwap(cipher, mask[i], 1 << (param.log2n - 1 - i), scheme, ring, bootHelper);	
        bootHelper.bootstrapping(cipher, param.logq, param.logQ, param.logT);
        timeutils.stop("CompAndSwap");
        
        CyclicArray ca(mvec, n);
        plainSort.compAndSwap(ca, mask, i, 1 << (param.log2n - 1 - i), true);
        mvec = ca.getArray();
        // for(int i = 0; i < n / 2; i++) {
        //     if(mvec[i] > mvec[i + n / 2]) {
        //         double x = mvec[i];
        //         mvec[i] = mvec[i + n / 2];
        //         mvec[i + n / 2] = x;
        //     }
        // }

        complex<double>* dvec = scheme.decrypt(secretKey, cipher);
        // PrintUtils::printArrays(mvec, dvec, n);
        PrintUtils::averageDifference(mvec, dvec, n);
    }
}

void TestBoot::testSelfBitonicMerge(Parameter param, long iter) {
	srand(time(NULL));
	SetNumThreads(8);
    
    long n = 1 << param.log2n;
	
    PrintUtils::parameter(param, "TestBoot::testSelfBitonicMerge");

    TimeUtils timeutils;
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


    double* mvec = EvaluatorUtils::randomRealArray(n);
	PlainSort ps;
	CyclicArray ca1(mvec, n);
	ps.runPlainSorting(ca1, param.log2n, false);
	Ciphertext cipher = scheme.encrypt(ca1.getArray(), n, param.logp, param.logQ);

    MaskingGenerator mg(param.log2n);
    double** mask = mg.getBitonicMergeMasking();

    timeutils.start("testSelfBitonicMerge");
    BootAlgo bootAlgo(param, iter);
    bootAlgo.selfBitonicMerge(cipher, mask, scheme, ring, bootHelper);
    timeutils.stop("testSelfBitonicMerge");

	CyclicArray ca(mvec, n);
	PlainSort plainSort;
	plainSort.selfBitonicMerge(ca, param.log2n, mask, true);

	complex<double>* dvec = scheme.decrypt(secretKey, cipher);
	PrintUtils::printArrays(ca.getArray(), dvec, n);
    PrintUtils::averageDifference(ca.getArray(), dvec, n);
}