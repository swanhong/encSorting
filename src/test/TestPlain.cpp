#include "TestPlain.h"

void TestPlain::plainSort(long logn) {
    CyclicArray ca;
    ca.randomGen(1 << logn);

    cout << "Before sorting" << endl;
    ca.printAsVector();
    
    cout << "After sorting" << endl;
    PlainSort plainSort;
    plainSort.runSorting(ca, logn);
    ca.printAsVector();
}

void TestPlain::showMasking(long log2n) {
    MaskingGenerator mg(log2n);
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
