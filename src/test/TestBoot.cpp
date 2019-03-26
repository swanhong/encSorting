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
	
	Ciphertext cipher = scheme.encrypt(mvec, n, parameter.logp, parameter.logQ);
	mvec = scheme.decrypt(secretKey, cipher);


	timeutils.start("Improved bootstrapping");
	boothelper.bootstrapping(cipher, parameter.logq, parameter.logQ, parameter.logT);
	timeutils.stop("Improved bootstrapping");

	cout << "* Before logQ = " << parameter.logq << endl;
	cout << "* After logQ = " << cipher.logq << endl;

	// Print Result and Difference //
    complex<double>* dvec = scheme.decrypt(secretKey, cipher);
	cout << "log2(avg of error) = " << diff(mvec, dvec, n) << endl;
	
	if(mvec != NULL) delete[] mvec;
	if(dvec != NULL) delete[] dvec;
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
    BootAlgo bootAlgo;
	bootAlgo.approxSqrt(cipher, parameter, iter, scheme, bootHelper);
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
    BootAlgo bootAlgo;
    bootAlgo.minMax(cipher1, cipher2, iter, param, scheme, boothelper);
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
    double** mask = mg.getMasking();

    timeutils.start("CompAndSwap");
    BootAlgo bootAlgo;
    bootAlgo.compAndSwap(cipher, mask[0], 1, iter, param, scheme, ring, bootHelper);
    timeutils.stop("CompAndSwap");

    for(int i = 0; i < n / 2; i++) {
        if(mvec[2 * i] > mvec[2 * i + 1]) {
            double x = mvec[2 * i];
            mvec[2 * i] = mvec[2 * i + 1];
            mvec[2 * i + 1] = x;
        }
    }
	complex<double>* dvec = scheme.decrypt(secretKey, cipher);
    PrintUtils::averageDifference(mvec, dvec, n);
}