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

void TestPlain::bitonicMerge(long log2n, long logNum) {
    PlainSort plainSort;
    CyclicArray* ca = new CyclicArray[1 << logNum];
    for(int i = 0; i < (1 << logNum); i++) {
        double* mvec = EvaluatorUtils::randomRealArray(1 << log2n);
        ca[i] = CyclicArray(mvec, 1 << log2n);
        plainSort.runPlainSorting(ca[i], log2n, (i % 2 == 0));
        ca[i].printAsVector();
        cout << "==" << endl;
    }

    plainSort.bitonicMerge(ca, log2n, logNum);
    for(int i = 0; i < (1 << logNum); i++) {
        ca[i].printAsVector();
        cout << "===" << endl;
    }
    
}

void TestPlain::plainTableSort(long logn, long logDataNum, bool increase) {
    CyclicArray ca;
    ca.randomGen(1 << logn);

    CyclicArray plain(ca);
    PlainSort plainSort;
    plainSort.runPlainTableSorting(ca, logn, logDataNum, 0);
    // ca.printAsVector();
    for (int i = 0; i < ca.length; i++) {
        cout << i << " : " << plain.get(i) << " // " << ca.get(i) << endl;
    }
    
}