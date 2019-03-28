#include "TestPlain.h"

void TestPlain::plainSort(long logn, bool increase) {
    CyclicArray ca;
    ca.randomGen(1 << logn);

    cout << "Before sorting" << endl;
    ca.printAsVector();
    
    cout << "After sorting" << endl;
    PlainSort plainSort;
    plainSort.runPlainSorting(ca, logn, increase);
    ca.printAsVector();
}

void TestPlain::showMasking(long log2n, bool increase) {
    MaskingGenerator mg(log2n, increase);
    double** mask = mg.getMasking();

    long length = 1 << log2n;
    
    for(int i = 0; i < (log2n + 1) * log2n / 2; i++){
        cout << "mask[" << i << "] = [";
        for(int j = 0; j < length; j++) {
            cout << mask[i][j] << ", ";
        }
        cout << "]" << endl;        
    }
}

void TestPlain::showBitonicMergeMasking(long log2n) {
    MaskingGenerator mg(log2n);
    double** mask = mg.getBitonicMergeMasking();

    long length = 1 << log2n;
    
    for(int i = 0; i < log2n; i++){
        cout << "mask[" << i << "] = [";
        for(int j = 0; j < length; j++) {
            cout << mask[i][j] << ", ";
        }
        cout << "]" << endl;        
    }
}
