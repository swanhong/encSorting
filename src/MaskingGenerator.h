#ifndef MASKINGGENERATOR_H_
#define MASKINGGENERATOR_H_

#include "iostream"

/*
 * class MaskingGenerator
 *      generate masking vectors for Batcher's Odd-Even Sort
 */
class MaskingGenerator {
private:
    double** mask;
    long log2n;
    long length;
    long logData;
    bool increase;
    

public:
    // The initializer automatically generates masking vectors
    MaskingGenerator(long _log2n, bool _increase = true, long _logData = 0);
    ~MaskingGenerator();

    long generateMaskingRec(long logNum, long logJump, long loc);

    void generateMaskingComparison(long loc, long jump);

    void generateMaskingMerge(long loc, long num, long jump);

    // outputs mask
    double** getMasking();

    double** getBitonicMergeMasking();
};

#endif // !MASKINGGENERATOR_H_