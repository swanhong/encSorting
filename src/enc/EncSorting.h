#ifndef ENCSORTING_H_
#define ENCSORTING_H_

#include "BootAlgorithm.h"
#include "EncAlgo.h"
#include "../plain/PlainSorting.h"
#include "../Parameter.h"
#include "../MaskingGenerator.h"

class EncSorting {
public:

    Parameter param;
    bool increase = true; 
    long inc;

    long logDataNum = 0;
    long colNum = 0;
    
    // Enc
    EncAlgo* encAlgo;
    ZZ*** maskPoly;
    bool maskPolyGen = false;
    ZZ*** maskRightPoly;
    ZZ*** maskTablePoly;
    ZZ*** maskTableRightPoly;
    bool maskTableGen = false;

    // Plain
    bool showDiff = false;

    PlainSort* plainSort;
    CyclicArray* ca;


    EncSorting(Parameter _param, long* iter, long numOfIter, bool = true, bool = false);
    ~EncSorting() {}

    Ciphertext encrypt(double* mvec);
    complex<double>* decrypt(Ciphertext& cipher);

    void showDiffFromPlain();

    void genMaskPoly();
    void genTableMaskPoly();

    void runEncSorting(Ciphertext& cipher);
    long sortingRecursion(Ciphertext& cipher, long logNum, long logDist, long loc);

    void runEncTableSorting(Ciphertext& cipher);

    long sortingTableRecursion(Ciphertext& cipher, long logNum, long logDist, long loc);

    void encSwapAndCompare(Ciphertext& cipher, long logDist, long loc);
    void encSwapTableAndCompare(Ciphertext& cipher, long dist, long loc);

    void bitonicMerge(Ciphertext* cipher, long logNum, BootScheme& scheme, Ring& ring, BootHelper& bootHelper);
    void bitonicMergeRec(Ciphertext* cipher, long start, long logNum, double** maskIncrease, double** maskDecrease, BootScheme& scheme, Ring& ring, BootHelper& bootHelper, bool increase);

    void bitonicTableMerge(Ciphertext* cipher, long logNum, long logDataNum, long colNum, BootScheme& scheme, Ring& ring, BootHelper& bootHelper, SecretKey& sk);
    void bitonicTableMergeRec(Ciphertext* cipher, long start, long logNum, long logDataNum, long colNum, ZZ** maskminMaxTablePoly, double*** maskInc, double*** maskDec, BootScheme& scheme, Ring& ring, BootHelper& bootHelper, SecretKey& sk, bool increase);
    
    void reverseHalf(Ciphertext* cipher, long logNum, BootScheme& scheme, Ring& ring, BootHelper& bootHelper);

    void runBitonicSort(Ciphertext& cipher, BootScheme& scheme, Ring& ring, BootHelper& bootHelper, bool increase);

    void runBitonicSortDec(Ciphertext& cipher, BootScheme& scheme, Ring& ring, BootHelper& bootHelper, bool increase, SecretKey sk);
};

#endif // !ENCSORTING_H_