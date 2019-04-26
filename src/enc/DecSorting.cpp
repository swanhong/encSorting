#include "DecSorting.h"

void DecSorting::runDecSorting(Parameter param, long iter, bool increase) {
    TimeUtils timeutils;
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

    //
    MaskingGenerator mg(param.log2n, increase);
    double** mask = mg.getMasking();
    
    bootAlgo = BootAlgo(param, iter, increase);
    //

    long n = 1 << param.log2n;
    double* mvec = EvaluatorUtils::randomRealArray(n);
    
    CyclicArray ca(mvec, n);
    
    Ciphertext cipher = scheme.encrypt(mvec, n, param.logp, param.logQ);

    decSortingRecursion(ca, cipher, param.log2n, 0, 0, mask, param, secretKey, scheme, ring, bootHelper);
}

long DecSorting::decSortingRecursion(CyclicArray& ca, Ciphertext& cipher, long logNum, long logJump, long loc, double** mask, Parameter& param, SecretKey& secretKey, BootScheme& scheme, Ring& ring, BootHelper& bootHelper) {
    TimeUtils timeutils;
    if (logNum == 1) {
        timeutils.start(to_string(loc)+"th CompAndSwap");
        compSwapPrintAll(ca, cipher, mask, loc, logJump, param, secretKey, scheme, ring, bootHelper);
        timeutils.stop(to_string(loc)+"th CompAndSwap");
    } else {
        if (logJump == 0) {
            loc = decSortingRecursion(ca, cipher, logNum - 1, logJump, loc, mask, param, secretKey, scheme, ring, bootHelper);
        }
        loc = decSortingRecursion(ca, cipher, logNum - 1, logJump + 1, loc, mask, param, secretKey, scheme, ring, bootHelper);

        timeutils.start(to_string(loc)+"th CompAndSwap");
        compSwapPrintAll(ca, cipher, mask, loc, logJump, param, secretKey, scheme, ring, bootHelper);
        timeutils.stop(to_string(loc)+"th CompAndSwap");
    }
    return loc + 1;
}

void DecSorting::compSwapPrintAll(CyclicArray& ca, Ciphertext& cipher, double** mask, long loc, long logJump, Parameter& param, SecretKey& secretKey, BootScheme& scheme, Ring& ring, BootHelper& bootHelper) {
    cout << "EncCompSwap" << endl;
    bootAlgo.compAndSwap(cipher, mask[loc], 1<< logJump, scheme, ring, bootHelper);
    cout << "PlainCompSwap" << endl;
    plainSort.compAndSwap(ca, mask, loc, 1 << logJump, true);
    cout << "print" << endl;
    decAndPrint(ca, cipher, param, scheme, secretKey);
}

void DecSorting::decAndPrint(CyclicArray& ca, Ciphertext& cipher, Parameter& param, BootScheme& scheme, SecretKey& secretKey) {
    complex<double>* dvec = scheme.decrypt(secretKey, cipher);
    double* mvec = ca.getArray();

    double max = 0;
    long max_index = 0;
    for (int i = 0; i < (1 << param.log2n); i++) {
        double diff = abs(mvec[i] - dvec[i].real());
        if (diff > max) {
            max = diff;
            max_index = i;
        }
        
        cout << i << " : " << mvec[i] << " // " << dvec[i].real() << " -- " << diff << endl;
    }
    // PrintUtils::printArrays(mvec, dvec, 1 << param.log2n);
    PrintUtils::averageDifference(mvec, dvec, 1<< param.log2n);
    cout << "log2(max error) = " << max_index << " : log2(" << max << ") = " << log2(max) << endl;
}