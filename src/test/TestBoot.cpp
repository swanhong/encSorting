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

void TestBoot::approxInverse(Parameter parameter, long iter) {
	long n = 1 << parameter.log2n;

	PrintUtils::parameter(parameter, "TestBoot::approxInverse");
	
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
    for (int i = 0; i < n; i++) {
        mvec[i] += 0.5;
    }
    
	Ciphertext cipher = scheme.encrypt(mvec, n, parameter.logp, parameter.logQ);

    
    timeutils.start("approxInverse");
    BootAlgo bootAlgo(parameter, iter);
	bootAlgo.approxInverse(cipher, scheme, bootHelper);
    timeutils.stop("approxInverse");

    // Print Result and Difference //	
	complex<double>* dvec = scheme.decrypt(secretKey, cipher);
    double* invvec = new double[n];
    for(int i = 0; i < n; i++) {
        invvec[i] = 1 / mvec[i];
    }
    for (int i = 0; i < n; i++)
    {
        cout << i << " : " << mvec[i] << " // " << invvec[i] << " // " << dvec[i].real() << endl;
    }
    
    PrintUtils::averageDifference(invvec, dvec, n);
}

void TestBoot::approxComp(Parameter parameter, long invIter, long compIter) {
    long n = 1 << parameter.log2n;

	PrintUtils::parameter(parameter, "TestBoot::approxComp");
	
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
    cout << 1 << endl;
	timeutils.stop("Bootstrapping Helper construct");


	double* mvec1 = EvaluatorUtils::randomRealArray(n);
    double* mvec2 = EvaluatorUtils::randomRealArray(n);
    for (int i = 0; i < n; i++) {     
        mvec1[i] += 0.5;
        mvec2[i] += 0.5;
                if (i < 10) {
        }
    }
    mvec1[52] = 0.7;

	Ciphertext cipher1 = scheme.encrypt(mvec1, n, parameter.logp, parameter.logQ);
    Ciphertext cipher2 = scheme.encrypt(mvec2, n, parameter.logp, parameter.logQ);

    
    timeutils.start("approxComp");
    BootAlgo bootAlgo(parameter, invIter, compIter);
	bootAlgo.comparison(cipher1, cipher2, scheme, bootHelper);
    timeutils.stop("approxComp");

    // Print Result and Difference //	
	complex<double>* dvec = scheme.decrypt(secretKey, cipher1);
    complex<double>* dvec2 = scheme.decrypt(secretKey, cipher2);
    double* compvec = new double[n];    
    for (int i = 0; i < n; i++) {
        if (mvec1[i] >= mvec2[i]) {
            compvec[i] = 1.;
        } else {
            compvec[i] = 0.;
        }
        cout << i << " : " << mvec1[i] << " // " << mvec2[i] << " // " << compvec[i] << " // "  << dvec[i].real() << " // " << dvec2[i].real() << endl;
    }
    
    PrintUtils::averageDifference(compvec, dvec, n);   
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

    // MaskingGenerator mg(param.log2n);
    // double** mask = mg.getMasking();

    // timeutils.start("CompAndSwap");
    // BootAlgo bootAlgo(param, iter);
    // bootAlgo.compAndSwap(cipher, mask[0], 1, scheme, ring, bootHelper);
	
    // timeutils.stop("CompAndSwap");

    // for(int i = 0; i < n / 2; i++) {
    //     if(mvec[2 * i] > mvec[2 * i + 1]) {
    //         double x = mvec[2 * i];
    //         mvec[2 * i] = mvec[2 * i + 1];
    //         mvec[2 * i + 1] = x;
    //     }
    // }

	MaskingGenerator mg(param.log2n);
    double** maskIncrease = mg.getBitonicMergeMasking();
	BootAlgo bootAlgo(param, iter);
	bootAlgo.compAndSwap(cipher, maskIncrease[0], 1 << (param.log2n - 1), scheme, ring, bootHelper);
	for(int i = 0; i < n / 2; i++) {
        if(mvec[i] > mvec[i + n / 2]) {
            double x = mvec[i];
            mvec[i] = mvec[i + n / 2];
            mvec[i + n / 2] = x;
        }
    }


	complex<double>* dvec = scheme.decrypt(secretKey, cipher);
	PrintUtils::printArrays(mvec, dvec, n);
    PrintUtils::averageDifference(mvec, dvec, n);
}

void TestBoot::reverse(Parameter param) {
    srand(time(NULL));
	SetNumThreads(8);
    
    long n = 1 << param.log2n;
	
    PrintUtils::parameter(param, "TestBoot::reverse");

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
	BootAlgo bootAlgo(param, 0);
	bootAlgo.reverse(cipher, mask, scheme, ring, bootHelper);

	complex<double>* dvec = scheme.decrypt(secretKey, cipher);
	PrintUtils::printArrays(mvec, dvec, n);
}

void TestBoot::compAndSwapTable(Parameter parameter, long logDataNum, long colNum, long invIter, long compIter) {
    long n = 1 << parameter.log2n;

	PrintUtils::parameter(parameter, "TestBoot::compAndSwapTable");
	
	TimeUtils timeutils;
	timeutils.start("KeyGen");
	Ring ring(parameter.logN, parameter.logQ);
	SecretKey secretKey(ring);
	BootScheme scheme(secretKey, ring);
	scheme.addConjKey(secretKey);
	scheme.addLeftRotKeys(secretKey);
    scheme.addRightRotKeys(secretKey);
	timeutils.stop("KeyGen");

	timeutils.start("Bootstrapping Helper construct");
	BootHelper bootHelper(parameter.log2n, parameter.radix, parameter.logc, scheme, ring, secretKey);
	timeutils.stop("Bootstrapping Helper construct");

	double* mvec = EvaluatorUtils::randomRealArray(n);
    // for (int i = 0; i < n; i++) {
    //     mvec[i] += 0.5;
    // }
    
	Ciphertext cipher = scheme.encrypt(mvec, n, parameter.logp, parameter.logQ);

    MaskingGenerator mg(parameter.log2n, logDataNum);
    double** mask = mg.getMasking();
    MaskingGenerator mg2(parameter.log2n, logDataNum);
    double** maskOther = mg2.getMaskingOther();
    MaskingGenerator mgTable(parameter.log2n, logDataNum, colNum, true);
    double** maskTable = mgTable.getMasking();
    MaskingGenerator mgTable2(parameter.log2n, logDataNum, colNum, true);
    double** maskTableOther = mgTable2.getMaskingOther();
    
    // long logn = parameter.log2n - logDataNum;
    // long maskNum = logn * (logn + 1) / 2;
    // PrintUtils::printSingleMatrix("mask", mask, maskNum, n);
    // cout << endl;
    // PrintUtils::printSingleMatrix("maskOther", maskOther, maskNum, n);
    // cout << endl;
    // PrintUtils::printSingleMatrix("maskTable", maskTable, maskNum, n);
    // cout << endl;
    // PrintUtils::printSingleMatrix("maskTableOther", maskTableOther, maskNum, n);
    
    timeutils.start("compAndSwapTable");
    BootAlgo bootAlgo(parameter, invIter, compIter);
	// bootAlgo.compAndSwapTable(cipher, logDataNum, mask[0], maskOther[0], maskTable[0], maskTableOther[0], 1 << logDataNum, scheme, ring, bootHelper);
    bootAlgo.compAndSwapTable(cipher, logDataNum, mask[1], maskOther[1], maskTable[1], maskTableOther[1], 1 << (logDataNum + 1), scheme, ring, bootHelper);
    timeutils.stop("compAndSwapTable");
    // Print Result and Difference //	
	complex<double>* dvec = scheme.decrypt(secretKey, cipher);
    // for (int i = 0; i < n; i++)
    // {
    //     cout << i << " : " << mvec[i] << " // " << invvec[i] << " // " << dvec[i].real() << endl;
    // }
    PrintUtils::printArrays(mvec, dvec, n);
    PrintUtils::averageDifference(mvec, dvec, n);
}