#ifndef DECSORTING_H_
#define DECSORTING_H_

#include "BootAlgorithm.h"
#include "EncSorting.h"
#include "../bootsrc/new_bootstrapping.h"
#include "../PrintUtils.h"
#include "../Parameter.h"
#include "EncAlgorithm.h"
#include "../plain/PlainSorting.h"

class DecSorting : public EncSorting {
    PlainSort plainSort;
public:
    void runDecSorting(Parameter param, long iter, bool increase);
    
    long decSortingRecursion(CyclicArray& ca, Ciphertext& cipher, long logNum, long logJump, long loc, double** mask, Parameter& param, SecretKey& secretKey, BootScheme& scheme, Ring& ring, BootHelper& bootHelper);
    
    void compSwapPrintAll(CyclicArray& ca, Ciphertext& cipher, double** mask, long loc, long logJump, Parameter& param, SecretKey& secretKey, BootScheme& scheme, Ring& ring, BootHelper& bootHelper);

    void decAndPrint(CyclicArray& ca, Ciphertext& cipher, Parameter& param, BootScheme& scheme, SecretKey& secretKey);
};

#endif //