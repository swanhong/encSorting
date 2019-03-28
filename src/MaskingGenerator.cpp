#include "MaskingGenerator.h"


MaskingGenerator::MaskingGenerator(long _log2n, bool _increase)  {
    log2n = _log2n;
    length = 1 << log2n;
    increase = _increase;
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
            if(increase) {
                mask[loc][i * jump * 2 + j] = 1;
            } else {
                mask[loc][length - 1 - (i * jump * 2 + j)] = 1;
            }
            
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
                if(increase) {
                    mask[loc][i * jump * num + (2 * k + 1) * jump + j] = 1;    
                } else {
                    mask[loc][length - 1 - (i * jump * num + (2 * k + 1) * jump + j)] = 1;    
                }
                
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
    long maskNum = log2n * (log2n + 1) / 2;
    mask = new double*[maskNum];
    for(int i = 0; i < maskNum; i++) {
        mask[i] = new double[length];
    }
    generateMaskingRec(log2n, 0, 0);
    return mask;
}

double** MaskingGenerator::getBitonicMergeMasking() {
    mask = new double*[log2n];
    for(int i = 0; i < log2n; i++) {
        mask[i] = new double[length];
        for(int j = 0; j < length; j++) {
            mask[i][j] = 0;
        }
        
    }

    for(int i = 0; i < log2n; i++) {
        for(int j = 0; j < (1 << i); j++) {
            for(int k = 0; k < (1 << (log2n - 1 - i)); k++) {
                mask[i][j * (1 << (log2n - i)) + k] = 1; 
            }
        }
    }
    return mask;
}
