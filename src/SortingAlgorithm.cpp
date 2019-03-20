#include "SortingAlgorithm.h"


SortingAlgorithm::SortingAlgorithm(CyclicArray ca_, long log2n_) {
    ca = ca_;
    log2n = log2n_; 
    maskNum = (log2n + 1) * log2n / 2;
    genAllMasking(log2n, mask);
}

void SortingAlgorithm::BatcherOddEvenSort() {
    BatcherOddEvenSortRec(log2n, 0, 0);
}

long SortingAlgorithm::BatcherOddEvenSortRec(long logNum, long logJump, long loc) {
    // cout << "call Rec(" << logNum << ", " << logJump << ")" << endl;
    if (logNum == 1) {
        // cout << "mask[" << loc << "] <- S(" << logNum << ", " << logJump << ")" << endl;
        // CompAndSwap()
        compAndSwap(loc, 1 << logJump);
        // mask[loc] = genMaskingComp(1 << log2n, 1 << logJump);
        return loc + 1;
    } else {
        if (logJump == 0) {
            loc = BatcherOddEvenSortRec(logNum - 1, logJump, loc);
        }
        loc = BatcherOddEvenSortRec(logNum - 1, logJump + 1, loc);
        // cout << "mask[" << loc << "] <- S(" << logNum << ", " << logJump << ")" << endl;
        // mask[loc] = genMaskingMerge(1 << log2n, 1 << logNum, 1 << logJump);
        compAndSwap(loc, 1 << logJump);
        return loc + 1;
    }
}

void SortingAlgorithm::compAndSwap(long loc, long dist) {
    long length = ca.length;
    CyclicArray maskCA(mask[loc], length);     
    
    CyclicArray dummy(length);
    mult(dummy, ca, maskCA);
    ca.sub(dummy);
    // cout << "ca" << endl;
    // ca.printAsVector();
    // cout << "maskCA" << endl;
    // maskCA.printAsVector();
    // cout << "dummy" << endl;
    // dummy.printAsVector();
    dummy.rightRotate(dist);
    
    getMinMax(dummy, ca);
    // std::cout << "after minmax" << '\n';
    // dummy.printAsVector();
    // ca.printAsVector();
    dummy.leftRotate(dist);
    ca.add(dummy);
}

double* genMaskingComp(long length, long jump) {   
    double* mask = new double[length];
    for(int i = 0; i < length; i++) {
        mask[i] = 0;
    }
    
    long repeat = length / (jump * 2);
    for(int i = 0; i < repeat; i++) {
        for(int j = 0; j < jump; j++) {
            mask[i * jump * 2 + j] = 1;
        }
    }
    return mask;
}

double* genMaskingMerge(long length, long num, long jump) {
    double* mask = new double[length];
    for(int i = 0; i < length; i++) {
        mask[i] = 0;
    }
    long repeat = length / (jump * num);
    for(int i = 0; i < repeat; i++) {
        for(int j = 0; j < jump; j++) {
            for(int k = 0; k < num / 2 - 1; k++) {
                mask[i * jump * num + (2 * k + 1) * jump + j] = 1;    
            }
            
            
        }        
    }
    return mask;
}


void genAllMasking(long log2n, double**& mask) {
    int maskNum = (log2n + 1) * log2n / 2;
    mask = new double*[maskNum];
    genMaskingRec(log2n, log2n, 0, mask, 0);
}

long genMaskingRec(long log2n, long logNum, long logJump, double**& mask, long loc) {

    // cout << "call Rec(" << log2n << ", " << logNum << ", " << logJump << ")" << endl;
    if (logNum == 1) {
        // cout << "mask[" << loc << "] <- S(" << logNum << ", " << logJump << ")" << endl;
        mask[loc] = genMaskingComp(1 << log2n, 1 << logJump);
        return loc + 1;
    } else {
        if (logJump == 0) {
            loc = genMaskingRec(log2n, logNum - 1, logJump, mask, loc);
        }
        loc = genMaskingRec(log2n, logNum - 1, logJump + 1, mask, loc);
        // cout << "mask[" << loc << "] <- S(" << logNum << ", " << logJump << ")" << endl;
        mask[loc] = genMaskingMerge(1 << log2n, 1 << logNum, 1 << logJump);
        return loc + 1;
    }
}

void CyclicArray::randomGen(long length_) {
    length = length_;
    start = 0;
    data = new double[length];
    for(int i = 0; i < length; i++) {
        data[i] = 0;
        data[i] = (double) rand() / RAND_MAX;
    }
}

