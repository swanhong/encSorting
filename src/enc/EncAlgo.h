#ifndef ENCALGO_H_
#define ENCALGO_H_

#include "EncInterface.h"

class EncAlgo : public EncInterface {
protected:
public:
    long* iter;
    
    long minMaxLoc = 0;
    long invLoc = 0;
    long compLoc = 1;

    long logDataNum = 0;
    long colNum = 0;
    bool table = false;

    EncAlgo(Parameter _param, long* _iter, long numOfIters, bool = false);

    void sqrtAlgorithm(Ciphertext& cipher);
    void minMaxAlgorithm(Ciphertext& minCipher, Ciphertext& maxCipher);
    void comparisonAlgorithm(Ciphertext& a, Ciphertext& b);

    void evalFcn(Ciphertext& cipher);    
    void approxSqrt(Ciphertext& cipher);
    void approxSqrt2(Ciphertext& cipher);

    void minMax(Ciphertext& minCipher, Ciphertext& maxCipher);
    void newMinMax(Ciphertext& minCipher, Ciphertext& maxCipher);
    void encSwap(Ciphertext& cipher, ZZ* mask, long dist, bool = true);

    void selfBitonicMerge(Ciphertext& cipher, ZZ** mask, bool = true);
    void reverse(Ciphertext& cipher, ZZ** mask);
    void reverse(Ciphertext& cipher, ZZ** maskLeft, ZZ** maskRight, long level, bool = true);
    void halfCleaner(Ciphertext& cipher, ZZ* mask, long dist, bool = true);

    void approxInverse(Ciphertext& cipher);
    void comparison(Ciphertext& cipher1, Ciphertext& cipher2);
    void newComparison(Ciphertext& cipher1, Ciphertext& cipher2);
    void minMaxTable(Ciphertext& minCipher, Ciphertext& maxCipher, Ciphertext& minCipherTable, Ciphertext& maxCipherTable, ZZ* mask, ZZ* maskTable);
    void encSwapTable(Ciphertext& cipher, ZZ* maskLeft, ZZ* maskRight, ZZ* maskTableLeft, ZZ* maskTableRight, long dist, bool = true);
};

#endif // !ENCALGO_H_
