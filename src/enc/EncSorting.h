#ifndef ENCSORTING_H_
#define ENCSORTING_H_

#include "BootAlgorithm.h"
#include "../Parameter.h"
#include "../MaskingGenerator.h"

class EncSorting {
public:
    EncSorting() {}
    ~EncSorting() {}
    
    void runSorting(Ciphertext& cipher, Parameter param, long iter, BootScheme& scheme, Ring& ring, BootHelper& bootHelper);

    long sortingRecursion(Ciphertext& cipher, long logNum, long logJump, long loc, Parameter param, long iter, double** mask, BootAlgo& bootAlgo, BootScheme& scheme, Ring& ring, BootHelper& bootHelper);

};

#endif // !ENCSORTING_H_