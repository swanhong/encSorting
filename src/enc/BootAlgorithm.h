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
    long sqrtIter;
    long invIter;
    long compIter;
    bool increase;
public:
    
    BootAlgo() {}
    BootAlgo(Parameter _param, long _sqrtiter, bool=true);
    BootAlgo(Parameter _param, long _inviter, long _compIter, bool=true);
    BootAlgo(Parameter _param, long _sqrtIter, long _invIter, long _compIter, bool=true);
    ~BootAlgo() {}

    void approxSqrt(Ciphertext& cipher, BootScheme& scheme, BootHelper& bootHelper);
    
    void approxInverse(Ciphertext& cipher, BootScheme& scheme, BootHelper& bootHelper); 

    void approxInverseWithDec(Ciphertext& cipher, BootScheme& scheme, BootHelper& bootHelper, SecretKey& secretKey);

    void minMax(Ciphertext& minCipher, Ciphertext& maxCipher, BootScheme& scheme, BootHelper& bootHelper);
    
    void comparison(Ciphertext& cipher1, Ciphertext& cipher2, BootScheme& scheme, BootHelper& bootHelper);

    void compAndSwap(Ciphertext& cipher, double* mask, long dist, BootScheme& scheme, Ring& ring, BootHelper& bootHelper);

    void selfBitonicMerge(Ciphertext& cipher, double** mask, BootScheme& scheme, Ring& ring, BootHelper& bootHelper);

    void reverse(Ciphertext& cipher, double** mask, BootScheme& scheme, Ring& ring, BootHelper& bootHelper);

    void compAndSwapTable(Ciphertext& cipher, long logDataNum, double* mask, double* maskOther, double* maskTable, double* maskTableOther, long dist, BootScheme& scheme, Ring& ring, BootHelper& bootHelper, SecretKey& secretKey);
};

#endif // !BootAlgorithm_H_