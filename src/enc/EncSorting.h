#ifndef ENCSORTING_H_
#define ENCSORTING_H_

#include "BootAlgorithm.h"
#include "../Parameter.h"
#include "../MaskingGenerator.h"

class EncSorting {
private:
    Parameter param;
    long iter;
    BootAlgo bootAlgo;

public:
    EncSorting(Parameter _param, long _iter) : param(_param), iter(_iter) {}
    ~EncSorting() {}
    
    void runEncSorting(Ciphertext& cipher, BootScheme& scheme, Ring& ring, BootHelper& bootHelper, bool increase);

    long sortingRecursion(Ciphertext& cipher, long logNum, long logJump, long loc, double** mask, BootScheme& scheme, Ring& ring, BootHelper& bootHelper);

    void bitonicMerge(Ciphertext* cipher, long logNum, BootScheme& scheme, Ring& ring, BootHelper& bootHelper);

    void bitonicMergeRec(Ciphertext* cipher, long start, long logNum, double** mask, double** mask2, BootScheme& scheme, Ring& ring, BootHelper& bootHelper, bool increase);
};

#endif // !ENCSORTING_H_