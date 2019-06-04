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
    
    // masking polynomials
    ZZ*** maskPoly;
    bool maskPolyGen = false;
    
    ZZ*** maskRightPoly;
    ZZ*** maskTablePoly;
    ZZ*** maskTableRightPoly;
    bool maskTableGen = false;

    ZZ*** maskMergePoly;
    ZZ*** maskMergeColPoly;
    ZZ*** maskMergeRightPoly;
    ZZ*** maskMergeRightColPoly;
    bool maskMergeGen = false;

    ZZ** maskMinMaxTablePoly;
    ZZ*** maskMinMaxTableRightPoly;
    bool maskMinMaxTableGen = false;

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
    void genMergeMaskingPoly();
    void genMinMaxTablePoly();
    void genTableMergeMaskingPoly();

    void setInc(bool _increase);

    void runEncSorting(Ciphertext& cipher);
    void runEncSorting(Ciphertext& cipher, bool _increase);
    long sortingRecursion(Ciphertext& cipher, long logNum, long logDist, long loc);

    void runEncTableSorting(Ciphertext& cipher);

    long sortingTableRecursion(Ciphertext& cipher, long logNum, long logDist, long loc);

    void encSwapAndCompare(Ciphertext& cipher, long logDist, long loc);
    void encSwapTableAndCompare(Ciphertext& cipher, long dist, long loc);

    void bitonicMerge(Ciphertext* cipher, long logNum);
    void bitonicMergeRec(Ciphertext* cipher, long start, long logNum, bool increase);

    void bitonicTableMerge(Ciphertext* cipher, long logNum);
    void bitonicTableMergeRec(Ciphertext* cipher, long start, long logNum, bool increase);
    
    void reverseHalf(Ciphertext* cipher, long logNum);

    void runBitonicSort(Ciphertext& cipher, BootScheme& scheme, Ring& ring, BootHelper& bootHelper, bool increase);

    void runBitonicSortDec(Ciphertext& cipher, BootScheme& scheme, Ring& ring, BootHelper& bootHelper, bool increase, SecretKey sk);
};

#endif // !ENCSORTING_H_