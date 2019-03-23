#ifndef BootAlgorithm_H_
#define BootAlgorithm_H_

#include "BootScheme.h"
#include "../bootsrc/new_bootstrapping.h"
#include "../PrintUtils.h"
#include "../Parameter.h"
#include "EncAlgorithm.h"


class BootAlgo {
public:
    BootAlgo() {}
    ~BootAlgo() {}

    void approxSqrt(Ciphertext& outCipher, const Ciphertext& inCipher, Parameter parameter, long iter, BootScheme& scheme, BootHelper& bootHelper);
    
    void minMax(Ciphertext& minCipher, Ciphertext& maxCipher, Ciphertext& input1, Ciphertext& input2, long iter, Parameter& param, BootScheme& scheme, BootHelper& bootHelper);

    void compAndSwap(Ciphertext& outCipher, const Ciphertext& inCipher, double* mask, long dist, long iter, Parameter& param, BootScheme& scheme, Ring& ring, BootHelper& bootHelper);

    
};

#endif // !BootAlgorithm_H_