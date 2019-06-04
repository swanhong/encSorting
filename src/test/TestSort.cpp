#include "TestSort.h"

void TestSort::sort(Parameter param, long _iter, bool increase) {
    long n = 1 << param.log2n;
    long iter[1] = {_iter};
    EncSorting encSorting(param, iter, 1, increase, true);
    encSorting.showDiffFromPlain();

    // double* mvec = EvaluatorUtils::randomRealArray(n);
    double* mvec = new double[n];
    for (int i = 0; i < n; i++) {
        mvec[i] = 1. / n * (double) i;
    }
    random_shuffle(&mvec[0], &mvec[n]);
       
	Ciphertext cipher = encSorting.encrypt(mvec);

    TimeUtils timeutils;
    timeutils.start("EncSort");
    encSorting.runEncSorting(cipher);
    timeutils.stop("EncSort"); 

    // *********************

    // run PlainSort
    CyclicArray ca(mvec, 1 << param.log2n);
    PlainSort plainSort(param.log2n);
    plainSort.runPlainSorting(ca, increase);
    mvec = ca.getArray();

    // Print Result and Difference //	
	complex<double>* dvec = encSorting.decrypt(cipher);
    PrintUtils::printArrays(mvec, dvec, n);
    PrintUtils::averageDifference(mvec, dvec, n);
}

void TestSort::tableSort(Parameter param, long logDataNum, long colNum, long invIter, long compIter, bool increase) {
    long n = 1 << param.log2n;
    long iter[4] = {logDataNum, colNum, invIter, compIter};
    EncSorting encSorting(param, iter, 4, increase);
    // encSorting.showDiffFromPlain();
    
    double* mvec = EvaluatorUtils::randomRealArray(n);    
	Ciphertext cipher = encSorting.encrypt(mvec);

    TimeUtils timeutils;
    timeutils.start("EncTableSort");
    encSorting.runEncTableSorting(cipher);
    timeutils.stop("EncTableSort"); 

    CyclicArray ca(mvec, 1 << param.log2n);
    PlainSort plainSort(param.log2n, logDataNum, colNum);
    plainSort.runPlainTableSorting(ca, increase);
    mvec = ca.getArray();
    
//     // Print Result and Difference //	
	complex<double>* dvec = encSorting.decrypt(cipher);
    PrintUtils::printArraysWithDataNum(mvec, dvec, n, logDataNum, colNum);
    PrintUtils::averageDifference(mvec, dvec, n);    
}

void TestSort:: merge(Parameter param, long _iter, long logNum) {
    long n = 1 << param.log2n;
    long iter[1] = {_iter};
    EncSorting encSorting(param, iter, 1);
    // encSorting.showDiffFromPlain();
    
    long num = 1 << logNum;
    double ** mvec = new double*[num];
    Ciphertext* cipher = new Ciphertext[num];
 
    // sort (as bitonic sequence)
    PlainSort plainSort(param.log2n);
    CyclicArray* ca = new CyclicArray[num];
    for(int i = 0; i < num; i++) {
        mvec[i] = EvaluatorUtils::randomRealArray(n);
        ca[i] = CyclicArray(mvec[i], n);
        plainSort.runPlainSorting(ca[i], i % 2 == 0);
        cipher[i] = encSorting.encrypt(ca[i].getArray());
    }

    TimeUtils timeutils;
    timeutils.start("Bitonic Merge");
    encSorting.bitonicMerge(cipher, logNum);
    timeutils.stop("Bitonic Merge"); 

    //* decrypt
    complex<double>** dvec = new complex<double>*[num];
    for(int i = 0; i < num; i++) {
         dvec[i] = encSorting.decrypt(cipher[i]);
    } 

    //* run plain bitonicMerge
    plainSort.bitonicMerge(ca, param.log2n, logNum);

    //* Print merged arrays
    for (int i = 0; i < num; i++) {
        PrintUtils::printArrays(ca[i].getArray(), dvec[i], n);
    }

    //* Print log2(avg of errors)
    for (int i = 0; i < num; i++) {
        PrintUtils::averageDifference(ca[i].getArray(), dvec[i], n);
    }    
}

void TestSort:: tableMerge(Parameter param, long logNum, long logDataNum, long colNum, long invIter, long compIter, bool increase) {
    long n = 1 << param.log2n;
    long num = 1 << logNum;
    long iter[4] = {logDataNum, colNum, invIter, compIter};
    EncSorting encSorting(param, iter, 4, increase);
    // encSorting.showDiffFromPlain();

    double ** mvec = new double*[num];
    Ciphertext* cipher = new Ciphertext[num];
    PlainSort plainSort(param.log2n, logDataNum, colNum);
    CyclicArray* ca = new CyclicArray[num];
    for(int i = 0; i < num; i++) {
        mvec[i] = EvaluatorUtils::randomRealArray(n); 
        ca[i] = CyclicArray(mvec[i], n);
        plainSort.runPlainTableSorting(ca[i], i % 2 == 0);
        cipher[i] = encSorting.encrypt(ca[i].getArray());
    }    
    TimeUtils timeutils;
    timeutils.start("Bitonic Table Merge");
    encSorting.bitonicTableMerge(cipher, logNum);
    timeutils.stop("Bitonic Table Merge"); 

    //* decrypt
    complex<double>** dvec = new complex<double>*[num];
    for(int i = 0; i < num; i++) {
         dvec[i] = encSorting.decrypt(cipher[i]);
    } 

    //* run plain bitonicMerge
    plainSort.bitonicTableMerge(ca, param.log2n, logNum, logDataNum, colNum);

    //* Print merged arrays
    for (int i = 0; i < num; i++) {
        PrintUtils::printArraysWithDataNum(ca[i].getArray(), dvec[i], n, logDataNum, colNum);
    }

    //* Print log2(avg of errors)
    for (int i = 0; i < num; i++) {
        PrintUtils::averageDifference(ca[i].getArray(), dvec[i], n);
    }    
}

void TestSort::sortAndMerge(Parameter param, long _iter, long logNum) {
    long n = 1 << param.log2n;
    long iter[1] = {_iter};
    EncSorting encSorting(param, iter, 1);
    // encSorting.showDiffFromPlain();
    
    long num = 1 << logNum;

    double** mvec = new double*[num];
    Ciphertext* cipher = new Ciphertext[num];
    for(int i = 0; i < num; i++) {
        mvec[i] = EvaluatorUtils::randomRealArray(n);
        cipher[i] = encSorting.encrypt(mvec[i]);
    }

    TimeUtils timeutils, timeTotal;
    
    timeTotal.start("Sort And Merge");

    for (int i = 0; i < num; i++) {
        timeutils.start("EncSort " + to_string(i));
        encSorting.runEncSorting(cipher[i], i % 2 == 0);
        timeutils.stop("EncSort " + to_string(i));
    }

    // timeutils.start("Reverse");
    // encSorting.reverseHalf(cipher, logNum, scheme, ring, bootHelper);
    // timeutils.stop("Reverse");

    timeutils.start("Bitonic Merge");
    encSorting.bitonicMerge(cipher, logNum);
    timeutils.stop("Bitonic Merge"); 

    timeTotal.stop("Total Sort And Merge");

    //* run plain bitonicMerge
    PlainSort plainSort(param.log2n);
    CyclicArray* ca = new CyclicArray[num];
    for (int i = 0; i < num; i++) {
        ca[i] = CyclicArray(mvec[i], n);
        plainSort.runPlainSorting(ca[i], i % 2 == 0);
    }
    plainSort.bitonicMerge(ca, param.log2n, logNum);

    //* decrypt
    complex<double>** dvec = new complex<double>*[num];
    for(int i = 0; i < num; i++) {
         dvec[i] = encSorting.decrypt(cipher[i]);
    } 

    //* Print full arrays
    for (int i = 0; i < num; i++) {
        PrintUtils::printArrays(ca[i].getArray(), dvec[i], n);
        cout << "===" << endl;
    }

    //* Print log2(avg of errors)
    for (int i = 0; i < num; i++) {
        PrintUtils::averageDifference(ca[i].getArray(), dvec[i], n);
    }    
}

// void TestSort::bitonicSort(Parameter param, long iter, bool increase) {
//     srand(time(NULL));
// 	SetNumThreads(16);
//     TimeUtils timeutils;
// 	PrintUtils::parameter(param, "TestBitonicSort");
//     long n = 1 << param.log2n;
//     long logp = param.logp;

//     timeutils.start("TestBitonicSort KeyGen");
//     Ring ring(param.logN, param.logQ);
//     SecretKey secretKey(ring);
//     BootScheme scheme(secretKey, ring);
//     scheme.addConjKey(secretKey);
//     scheme.addLeftRotKeys(secretKey);
//     scheme.addRightRotKeys(secretKey);
//     timeutils.stop("TestBitonicSort KeyGen");

// 	timeutils.start("Bootstrapping Helper construct");
// 	BootHelper bootHelper(param.log2n, param.radix, param.logc, scheme, ring, secretKey);
// 	timeutils.stop("Bootstrapping Helper construct");

//     // double* mvec = EvaluatorUtils::randomRealArray(n);
//     double* mvec = new double[n];
//     for (int i = 0; i < n; i++) {
//         mvec[i] = 1. / n * (double) i;
//     }
//     random_shuffle(&mvec[0], &mvec[n]);
       
// 	Ciphertext cipher = scheme.encrypt(mvec, n, param.logp, param.logQ);

//     timeutils.start("EncSort");
//     EncSorting encSorting(param, iter);
//     encSorting.runBitonicSortDec(cipher, scheme, ring, bootHelper, increase, secretKey);
//     timeutils.stop("EncSort"); 

//     // run PlainSort
//     CyclicArray ca(mvec, 1 << param.log2n);
//     PlainSort plainSort;
//     plainSort.runPlainSorting(ca, param.log2n, increase);
//     mvec = ca.getArray();

//     // Print Result and Difference //	
// 	complex<double>* dvec = scheme.decrypt(secretKey, cipher);
//     PrintUtils::printArrays(mvec, dvec, n);
//     PrintUtils::averageDifference(mvec, dvec, n);
// }