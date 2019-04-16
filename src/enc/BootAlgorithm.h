#ifndef BootAlgorithm_H_
#define BootAlgorithm_H_

#include "BootScheme.h"
#include "../bootsrc/new_bootstrapping.h"
#include "../PrintUtils.h"
#include "../Parameter.h"
#include "EncAlgorithm.h"


class BootAlgo {
private:
    Parameter param;
    long iter;
    bool increase;
public:
    
    BootAlgo() {}
    BootAlgo(Parameter _param, long _iter, bool=true);
    ~BootAlgo() {}

    void approxSqrt(Ciphertext& cipher, BootScheme& scheme, BootHelper& bootHelper);
    
    void minMax(Ciphertext& minCipher, Ciphertext& maxCipher, BootScheme& scheme, BootHelper& bootHelper);

    void compAndSwap(Ciphertext& cipher, double* mask, long dist, BootScheme& scheme, Ring& ring, BootHelper& bootHelper);

    void selfBitonicMerge(Ciphertext& cipher, double** mask, BootScheme& scheme, Ring& ring, BootHelper& bootHelper);
};

#endif // !BootAlgorithm_H_