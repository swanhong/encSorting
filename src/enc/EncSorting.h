#ifndef ENCSORTING_H_
#define ENCSORTING_H_

#include "BootAlgorithm.h"
#include "EncAlgo.h"
#include "../plain/PlainSorting.h"
#include "../Parameter.h"
#include "../MaskingGenerator.h"

class EncSorting {
protected:
    Parameter param;
    EncAlgo* encAlgo;
    bool increase = true; 

    ZZ** maskSorting;
    bool maskSortingGen = false;
    ZZ** maskSortingDec;
    bool maskSortingDecGen = false;

public:
    EncSorting(Parameter _param, long* iter, long numOfIter, bool = true, bool = false);
    ~EncSorting() {}

    Ciphertext encrypt(double* mvec);
    complex<double>* decrypt(Ciphertext& cipher);

    void genMaskSorting();
    void genMaskSortingDec();
    
    void runEncSorting(Ciphertext& cipher, bool=true);
    long sortingRecursion(Ciphertext& cipher, long logNum, long logJump, long loc);

    void runEncTableSorting(Ciphertext& cipher, long logDataNum, long colNum, BootScheme& scheme, Ring& ring, BootHelper& bootHelper, SecretKey& secretKey, bool=true);
    
    long sortingRecursion(Ciphertext& cipher, long logNum, long logJump, long loc, double** mask, BootScheme& scheme, Ring& ring, BootHelper& bootHelper, SecretKey sk, PlainSort ps);
    
    long sortingTableRecursion(Ciphertext& cipher, long logDataNum, long colNum, long logNum, long logJump, long loc,
                                    double** mask, double** maskRight, double** maskTable, double** maskTableRight,
                                    BootScheme& scheme, Ring& ring, BootHelper& bootHelper, SecretKey& secretKey, PlainSort ps);

    void compAndSwapBothWithDec(Ciphertext& cipher, long logJump, long loc, double** mask, BootScheme& scheme, Ring& ring, BootHelper& bootHelper, SecretKey sk, PlainSort ps);
    void CompAndSwapTableBothWithDec(Ciphertext& cipher, long logDataNum, long colNum, long dist,
                                    double* mask, double* maskRight, double* maskTable, double* maskTableRight,
                                    BootScheme& scheme, Ring& ring, BootHelper& bootHelper, SecretKey& secretKey, PlainSort ps);

    void bitonicMerge(Ciphertext* cipher, long logNum, BootScheme& scheme, Ring& ring, BootHelper& bootHelper);
    void bitonicMergeRec(Ciphertext* cipher, long start, long logNum, double** maskIncrease, double** maskDecrease, BootScheme& scheme, Ring& ring, BootHelper& bootHelper, bool increase);

    void bitonicTableMerge(Ciphertext* cipher, long logNum, long logDataNum, long colNum, BootScheme& scheme, Ring& ring, BootHelper& bootHelper, SecretKey& sk);
    void bitonicTableMergeRec(Ciphertext* cipher, long start, long logNum, long logDataNum, long colNum, ZZ** maskminMaxTablePoly, double*** maskInc, double*** maskDec, BootScheme& scheme, Ring& ring, BootHelper& bootHelper, SecretKey& sk, bool increase);
    
    void reverseHalf(Ciphertext* cipher, long logNum, BootScheme& scheme, Ring& ring, BootHelper& bootHelper);

    void runBitonicSort(Ciphertext& cipher, BootScheme& scheme, Ring& ring, BootHelper& bootHelper, bool increase);

    void runBitonicSortDec(Ciphertext& cipher, BootScheme& scheme, Ring& ring, BootHelper& bootHelper, bool increase, SecretKey sk);
};

#endif // !ENCSORTING_H_