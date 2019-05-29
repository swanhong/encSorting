#include "TestAlgo.h"

void TestAlgo::bootstrapping(Parameter param) {
	cout << "Run TestAlgo::bootstrapping" << endl;
    
    long n = 1 << param.log2n;
    EncAlgo encAlgo(param, 0, 0);
    double* mvec = EvaluatorUtils::randomRealArray(n);
    Ciphertext cipher = encAlgo.encrypt(mvec);

    encAlgo.bootstrapping(cipher);

    complex<double>* dvec = encAlgo.decrypt(cipher);
	PrintUtils::printArrays(mvec, dvec, n);
    PrintUtils::averageDifference(mvec, dvec, n);
}

void TestAlgo::approxSqrt(Parameter parameter, long _iter) {
    cout << "Run TestAlgo::approxSqrt" << endl;

	long n = 1 << parameter.log2n;

    long iter[1] = {_iter};
    EncAlgo encAlgo(parameter, iter, 1);

	double* mvec = EvaluatorUtils::randomRealArray(n);
    
    Ciphertext cipher = encAlgo.encrypt(mvec);
    
    encAlgo.sqrtAlgorithm(cipher);

    complex<double>* dvec2 = encAlgo.decrypt(cipher);
    for(int i = 0; i < n; i++) {
        mvec[i] = sqrt(mvec[i]);
    }

    // PrintUtils::printArrays(mvec, dvec2, n);
    PrintUtils::averageDifference(mvec, dvec2, n);
}

void TestAlgo::minMax(Parameter param, long _iter) {
    cout << "Run TestAlgo::minMax" << endl;
    
    long n = 1 << param.log2n;

    long iter[1] = {_iter};
    EncAlgo encAlgo(param, iter, 1);

	double* mvec1 = EvaluatorUtils::randomRealArray(n);
    double* mvec2 = EvaluatorUtils::randomRealArray(n);
	Ciphertext cipher1 = encAlgo.encrypt(mvec1);
    Ciphertext cipher2 = encAlgo.encrypt(mvec2);
    
    TimeUtils timeutils;
    timeutils.start("minMax");
    // encAlgo.minMax(cipher1, cipher2);
    encAlgo.newMinMax(cipher1, cipher2);
    timeutils.stop("minMax");
    
    complex<double>* dmax = encAlgo.decrypt(cipher2);
    complex<double>* dmin = encAlgo.decrypt(cipher1);
    
    double* max = new double[n];
    double* min = new double[n];
    for(int i = 0; i < n; i++) {
        if(mvec1[i] > mvec2[i]) {
            max[i] = mvec1[i];
            min[i] = mvec2[i];
        } else {
            max[i] = mvec2[i];
            min[i] = mvec1[i];
        }   
    }
    // PrintUtils::printArrays(max, dmax, n);
    // PrintUtils::printArrays(min, dmin, n);
    PrintUtils::averageDifference(max, dmax, n);
    PrintUtils::averageDifference(min, dmin, n);
}

void TestAlgo::EncSwap(Parameter param, long _iter) {
    cout << "Run TestAlgo::EncSwap" << endl;
    
    long n = 1 << param.log2n;

    long iter[1] = {_iter};
    EncAlgo encAlgo(param, iter, 1);

    double* mvec = EvaluatorUtils::randomRealArray(n);

    Ciphertext cipher = encAlgo.encrypt(mvec);

    MaskingGenerator mg(param.log2n);
    double** mask = mg.getMasking();
    ZZ* maskPoly = encAlgo.encode(mask[0]);

    PrintUtils::printSingleArray("mask[0]", mask[0], n);

    TimeUtils timeutils;
    timeutils.start("EncSwap");
    encAlgo.encSwap(cipher, maskPoly, 1);
    timeutils.stop("EncSwap");

    for(int i = 0; i < n / 2; i++) {
        if(mvec[2 * i] > mvec[2 * i + 1]) {
            double x = mvec[2 * i];
            mvec[2 * i] = mvec[2 * i + 1];
            mvec[2 * i + 1] = x;
        }
    }
	complex<double>* dvec = encAlgo.decrypt(cipher);
	PrintUtils::printArrays(mvec, dvec, n);
    PrintUtils::averageDifference(mvec, dvec, n);
}

void TestAlgo::reverse(Parameter param) {
    cout << "Run TestAlgo::reverse" << endl;
    
    long n = 1 << param.log2n;

    EncAlgo encAlgo(param, 0, 0);

    double* mvec = EvaluatorUtils::randomRealArray(n);

    Ciphertext cipher = encAlgo.encrypt(mvec);

	MaskingGenerator mg(param.log2n);
    double** mask = mg.getBitonicMergeMasking();
    ZZ** maskPoly = new ZZ*[param.log2n];
    for (int i = 0; i < param.log2n; i++) {
        maskPoly[i] = encAlgo.encode(mask[i]);
    }
	
    encAlgo.reverse(cipher, maskPoly);
	// bootAlgo.reverse(cipher, mask, scheme, ring, bootHelper);

    complex<double>* dvec = encAlgo.decrypt(cipher);
	PrintUtils::printArrays(mvec, dvec, n);
}

void TestAlgo::halfCleaner(Parameter param, long iter) {
    cout << "Run TestAlgo::halfCleaner" << endl;
    
    long n = 1 << param.log2n;

    EncAlgo encAlgo(param, 0, 0);

    double* mvec = EvaluatorUtils::randomRealArray(n);

    Ciphertext cipher = encAlgo.encrypt(mvec);
	MaskingGenerator mg(param.log2n);
    double** mask = mg.getBitonicMergeMasking();

    ZZ* maskPoly = encAlgo.encode(mask[0]);

    PrintUtils::printSingleArray("mask[0]", mask[0], n);

    encAlgo.halfCleaner(cipher, maskPoly, 1 << (param.log2n - 1));

	complex<double>* dvec = encAlgo.decrypt(cipher);
	PrintUtils::printArrays(mvec, dvec, n);
}

void TestAlgo::approxInverse(Parameter param, long _iter) {
	cout << "Run TestTableAlgo::approxInverse" << endl;
    
    long n = 1 << param.log2n;
    long iter[1] = {_iter};
    EncAlgo encAlgo(param, iter, 1);
	double* mvec = EvaluatorUtils::randomRealArray(n);
    for (int i = 0; i < n; i++) {
        mvec[i] += 0.5;
        cout << "mvec[i]" << mvec[i] << endl;
    }

    Ciphertext cipher = encAlgo.encrypt(mvec);

    TimeUtils timeutils;   
    timeutils.start("approxInverse");
    long start = cipher.logq;
    encAlgo.approxInverse(cipher);
    timeutils.stop("approxInverse");
    cout << "consumed logq = " << start - cipher.logq << endl;
    
    // Print Result and Difference //	
	complex<double>* dvec = encAlgo.decrypt(cipher);
    double* invvec = new double[n];
    for(int i = 0; i < n; i++) {
        invvec[i] = 1 / mvec[i];
    }
    cout << "mvec // invvec // dvec" << endl;
    for (int i = 0; i < n; i++)
    {
        cout << i << " : " << mvec[i] << " // " << invvec[i] << " // " << dvec[i].real() << endl;
    }
    
    PrintUtils::averageDifference(invvec, dvec, n);
}

void TestAlgo::comparison(Parameter param, long invIter, long compIter) {
    cout << "Run TestAlgo::comparison" << endl;
    
    long n = 1 << param.log2n;

    long iter[2] = {invIter, compIter};
    EncAlgo encAlgo(param, iter, 2);

	double* mvec1 = EvaluatorUtils::randomRealArray(n);
    double* mvec2 = EvaluatorUtils::randomRealArray(n);
	Ciphertext cipher1 = encAlgo.encrypt(mvec1);
    Ciphertext cipher2 = encAlgo.encrypt(mvec2);
    
    TimeUtils timeutils;
    timeutils.start("approxComp");
	encAlgo.comparison(cipher1, cipher2);
    timeutils.stop("approxComp");

    // Print Result and Difference //	
	complex<double>* dvec1 = encAlgo.decrypt(cipher1);
	complex<double>* dvec2 = encAlgo.decrypt(cipher2);
    
    double* compvec = new double[n];    
    cout << "mvec1 // mvec2 // compvec // dvec1 // dvec2" << endl;
    for (int i = 0; i < n; i++) {
        if (mvec1[i] >= mvec2[i]) {
            compvec[i] = 1.;
        } else {
            compvec[i] = 0.;
        }
        cout << i << " : " << mvec1[i] << " // " << mvec2[i] << " // " << compvec[i] << " // "  << dvec1[i].real() << " // " << dvec2[i].real() << endl;
    }
    
    PrintUtils::averageDifference(compvec, dvec1, n);   
}

void TestAlgo::encSwapTable(Parameter param, long logDataNum, long colNum, long invIter, long compIter) {
    cout << "start TestAlgo::encSwapTable" << endl;
    long n = 1 << param.log2n;

    long iter[4] = {logDataNum, colNum, invIter, compIter};
	EncAlgo encAlgo(param, iter, 4);

	double* mvec = EvaluatorUtils::randomRealArray(n);
	Ciphertext cipher = encAlgo.encrypt(mvec);

    MaskingGenerator mg(param.log2n, logDataNum);
    double** mask = mg.getMasking();
    MaskingGenerator mg2(param.log2n, logDataNum);
    double** maskOther = mg2.getMaskingOther();
    MaskingGenerator mgTable(param.log2n, logDataNum, colNum, true);
    double** maskTable = mgTable.getMasking();
    MaskingGenerator mgTable2(param.log2n, logDataNum, colNum, true);
    double** maskTableOther = mgTable2.getMaskingOther();

    ZZ* maskLeft = encAlgo.encode(mask[0]);
    ZZ* maskRight = encAlgo.encode(maskOther[0]);
    ZZ* maskTableLeft = encAlgo.encode(maskTable[0]);
    ZZ* maskTableRight = encAlgo.encode(maskTableOther[0]);

    TimeUtils timeutils;
    timeutils.start("encSwapTable");
    encAlgo.encSwapTable(cipher, maskLeft, maskRight, maskTableLeft, maskTableRight, 1 << logDataNum);
    timeutils.stop("encSwapTable");

    PlainSort plainSort;
    CyclicArray ca(mvec, n);
    plainSort.compAndSwapTable(ca, logDataNum, colNum, mask[0], maskOther[0], maskTable[0], maskTableOther[0], 1 << (logDataNum), true);
    mvec = ca.getArray();
    
    // Print Result and Difference //	
	complex<double>* dvec = encAlgo.decrypt(cipher);
    PrintUtils::printArrays(mvec, dvec, n);
    PrintUtils::averageDifference(mvec, dvec, n);
}