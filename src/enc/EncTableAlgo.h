// #ifndef ENCTABLEALGO_H_
// #define ENCTABLEALGO_H_

// #include "../../HEAAN/src/HEAAN.h"
// #include "BootScheme.h"
// #include "EncInterface.h"
// #include "../Parameter.h"
// #include "../PrintUtils.h"


// class EncTableAlgo : public EncInterface {
// protected:
// public:
//     Parameter param;
//     Ring* ring;
//     SecretKey* secretKey;
//     BootScheme* scheme;
//     BootHelper* bootHelper;

//     long invIter;
//     long compIter;

//     long logDataNum;
//     long colNum;

//     bool printCondition;
    
//     EncTableAlgo(Parameter _param, long _logDataNum, long _colNum, long _invIter, long _compIter, bool = false);

//     void approxInverse(Ciphertext& cipher);

//     void comparison(Ciphertext& a, Ciphertext& b);

//     void minMaxTable(Ciphertext& minCipher, Ciphertext& maxCipher, Ciphertext& minCipherTable, Ciphertext& maxCipherTable, ZZ* mask, ZZ* maskTable);

//     void compAndSwapTable(Ciphertext& cipher, long logDataNum, long colNum, ZZ* maskLeft, ZZ* maskRight, ZZ* maskTableLeft, ZZ* maskTableRight, long dist, bool = true);

    
// };


// #endif //!ENCTABLEALGO_H_