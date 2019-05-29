#include "TestTableAlgo.h"

void TestTableAlgo::approxInverse(Parameter param, long iter) {
	cout << "Run TestTableAlgo::approxInverse" << endl;
    
    long n = 1 << param.log2n;
    EncTableAlgo encTableAlgo(param, 0, 0, iter, 0);
	double* mvec = EvaluatorUtils::randomRealArray(n);
    for (int i = 0; i < n; i++) {
        mvec[i] += 0.5;
        cout << "mvec[i]" << mvec[i] << endl;
    }

    Ciphertext cipher = encTableAlgo.encrypt(mvec);

    TimeUtils timeutils;   
    timeutils.start("approxInverse");
    long start = cipher.logq;
    encTableAlgo.approxInverse(cipher);
    timeutils.stop("approxInverse");
    cout << "consumed logq = " << start - cipher.logq << endl;
    
    // Print Result and Difference //	
	complex<double>* dvec = encTableAlgo.decrypt(cipher);
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

// void TestTableAlgo::approxComp(Parameter parameter, long invIter, long compIter) {
//     long n = 1 << parameter.log2n;

// 	PrintUtils::parameter(parameter, "TestAlgo::approxComp");
	
// 	TimeUtils timeutils;
// 	timeutils.start("KeyGen");
// 	Ring ring(parameter.logN, parameter.logQ);
// 	SecretKey secretKey(ring);
// 	BootScheme scheme(secretKey, ring);
// 	scheme.addConjKey(secretKey);
// 	scheme.addLeftRotKeys(secretKey);
// 	timeutils.stop("KeyGen");

// 	timeutils.start("Bootstrapping Helper construct");
// 	BootHelper bootHelper(parameter.log2n, parameter.radix, parameter.logc, scheme, ring, secretKey);
// 	timeutils.stop("Bootstrapping Helper construct");


// 	double* mvec1 = EvaluatorUtils::randomRealArray(n);
//     double* mvec2 = EvaluatorUtils::randomRealArray(n);

// 	Ciphertext cipher1 = scheme.encrypt(mvec1, n, parameter.logp, parameter.logQ);
//     Ciphertext cipher2 = scheme.encrypt(mvec2, n, parameter.logp, parameter.logQ);

    
//     timeutils.start("approxComp");
//     BootAlgo bootAlgo(parameter, invIter, compIter);
//     // bootHelper.bootstrapping(cipher1, parameter.logq, parameter.logQ, parameter.logT);
//     // bootHelper.bootstrapping(cipher2, parameter.logq, parameter.logQ, parameter.logT);
// 	bootAlgo.comparison(cipher1, cipher2, scheme, bootHelper);
//     timeutils.stop("approxComp");

//     // Print Result and Difference //	
// 	complex<double>* dvec = scheme.decrypt(secretKey, cipher1);
//     complex<double>* dvec2 = scheme.decrypt(secretKey, cipher2);
//     double* compvec = new double[n];    
//     cout << "mvec1 // mvec2 // compvec // dvec1 // dvec2" << endl;
//     for (int i = 0; i < n; i++) {
//         if (mvec1[i] >= mvec2[i]) {
//             compvec[i] = 1.;
//         } else {
//             compvec[i] = 0.;
//         }
//         cout << i << " : " << mvec1[i] << " // " << mvec2[i] << " // " << compvec[i] << " // "  << dvec[i].real() << " // " << dvec2[i].real() << endl;
//     }
    
//     PrintUtils::averageDifference(compvec, dvec, n);   
// }

// void TestTablelgo::compAndSwapTable(Parameter parameter, long logDataNum, long colNum, long invIter, long compIter) {
//     long n = 1 << parameter.log2n;

// 	PrintUtils::parameter(parameter, "TestAlgo::compAndSwapTable");
	
// 	TimeUtils timeutils;
// 	timeutils.start("KeyGen");
// 	Ring ring(parameter.logN, parameter.logQ);
// 	SecretKey secretKey(ring);
// 	BootScheme scheme(secretKey, ring);
// 	scheme.addConjKey(secretKey);
// 	scheme.addLeftRotKeys(secretKey);
//     scheme.addRightRotKeys(secretKey);
// 	timeutils.stop("KeyGen");

// 	timeutils.start("Bootstrapping Helper construct");
// 	BootHelper bootHelper(parameter.log2n, parameter.radix, parameter.logc, scheme, ring, secretKey);
// 	timeutils.stop("Bootstrapping Helper construct");

// 	double* mvec = EvaluatorUtils::randomRealArray(n);
//     // for (int i = 0; i < n; i++) {
//     //     mvec[i] += 0.5;
//     // }
    
// 	Ciphertext cipher = scheme.encrypt(mvec, n, parameter.logp, parameter.logQ);

//     MaskingGenerator mg(parameter.log2n, logDataNum);
//     double** mask = mg.getMasking();
//     MaskingGenerator mg2(parameter.log2n, logDataNum);
//     double** maskOther = mg2.getMaskingOther();
//     MaskingGenerator mgTable(parameter.log2n, logDataNum, colNum, true);
//     double** maskTable = mgTable.getMasking();
//     MaskingGenerator mgTable2(parameter.log2n, logDataNum, colNum, true);
//     double** maskTableOther = mgTable2.getMaskingOther();
    
//     long logn = parameter.log2n - logDataNum;
//     long maskNum = logn * (logn + 1) / 2;
//     // PrintUtils::printSingleMatrix("mask", mask, maskNum, n);
//     // cout << endl;
//     // PrintUtils::printSingleMatrix("maskOther", maskOther, maskNum, n);
//     // cout << endl;
//     // PrintUtils::printSingleMatrix("maskTable", maskTable, maskNum, n);
//     // cout << endl;
//     // PrintUtils::printSingleMatrix("maskTableOther", maskTableOther, maskNum, n);
    
//     complex<double>* dvec;

//     timeutils.start("compAndSwapTable");
//     BootAlgo bootAlgo(parameter, invIter, compIter);
//     // bootHelper.bootstrapping(cipher, parameter.logq, parameter.logQ, parameter.logT);
// 	bootAlgo.compAndSwapTable(cipher, logDataNum, colNum, mask[0], maskOther[0], maskTable[0], maskTableOther[0], 1 << logDataNum, scheme, ring, bootHelper, secretKey);
//     dvec = scheme.decrypt(secretKey, cipher);
//     // PrintUtils::printArrays(mvec, dvec, n);
//     PrintUtils::averageDifference(mvec, dvec, n);
//     bootAlgo.compAndSwapTable(cipher, logDataNum, colNum, mask[1], maskOther[1], maskTable[1], maskTableOther[1], 1 << (logDataNum + 1), scheme, ring, bootHelper, secretKey);
//     dvec = scheme.decrypt(secretKey, cipher);
//     // PrintUtils::printArrays(mvec, dvec, n);
//     PrintUtils::averageDifference(mvec, dvec, n);
//     bootAlgo.compAndSwapTable(cipher, logDataNum, colNum, mask[2], maskOther[2], maskTable[2], maskTableOther[2], 1 << (logDataNum), scheme, ring, bootHelper, secretKey);
//     dvec = scheme.decrypt(secretKey, cipher);
//     // PrintUtils::printArrays(mvec, dvec, n);
//     PrintUtils::averageDifference(mvec, dvec, n);
//     timeutils.stop("compAndSwapTable");
//     // Print Result and Difference //	
// 	// dvec = scheme.decrypt(secretKey, cipher);
//     // PrintUtils::printArrays(mvec, dvec, n);
//     // PrintUtils::averageDifference(mvec, dvec, n);
// }
