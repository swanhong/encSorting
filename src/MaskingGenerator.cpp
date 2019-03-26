#include "MaskingGenerator.h"


MaskingGenerator::MaskingGenerator(long _log2n)  {
    log2n = _log2n;
    length = 1 << log2n;

    long maskNum = log2n * (log2n + 1) / 2;
    mask = new double*[maskNum];
    for(int i = 0; i < maskNum; i++) {
        mask[i] = new double[length];
    }

    generateMaskingRec(log2n, 0, 0);
}

MaskingGenerator::~MaskingGenerator() {
    if(mask != NULL) delete[] mask;
}

void MaskingGenerator::generateMaskingComparison(long loc, long jump) {
    for(int i = 0; i < length; i++) {
        mask[loc][i] = 0;
    }
    
    long repeat = length / (jump * 2);
    for(int i = 0; i < repeat; i++) {
        for(int j = 0; j < jump; j++) {
            mask[loc][i * jump * 2 + j] = 1;
        }
    }
}

void MaskingGenerator::generateMaskingMerge(long loc, long num, long jump) {
    for(int i = 0; i < length; i++) {
        mask[loc][i] = 0;
    }
    long repeat = length / (jump * num);
    for(int i = 0; i < repeat; i++) {
        for(int j = 0; j < jump; j++) {
            for(int k = 0; k < num / 2 - 1; k++) {
                mask[loc][i * jump * num + (2 * k + 1) * jump + j] = 1;    
            }            
        }        
    }
}

long MaskingGenerator::generateMaskingRec(long logNum, long logJump, long loc) {
    if (logNum == 1) {
        generateMaskingComparison(loc, 1 << logJump);
    } else {
        if (logJump == 0) {
            loc = generateMaskingRec(logNum - 1, logJump, loc);
        }
        loc = generateMaskingRec(logNum - 1, logJump + 1, loc);
        generateMaskingMerge(loc, 1 << logNum, 1 << logJump);
    }
    return loc + 1;
}

double** MaskingGenerator::getMasking() {
    return mask;
}
