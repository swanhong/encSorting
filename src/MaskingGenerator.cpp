#include "MaskingGenerator.h"


MaskingGenerator::MaskingGenerator(long _log2n, bool _increase)  {
    log2n = _log2n;
    length = 1 << log2n;
    increase = _increase;
    logDataNum = 0;
    dataNum = 1;
    colNum = -1;
    setMaskNum();
}

MaskingGenerator::MaskingGenerator(long _log2n, long _logDataNum, bool _increase) {
    log2n = _log2n;
    length = 1 << log2n;
    increase = _increase;
    logDataNum = _logDataNum;
    dataNum = 1 << logDataNum;
    colNum = -1;
    setMaskNum();
}

MaskingGenerator::MaskingGenerator(long _log2n, long _logDataNum, long _colNum, bool _increase) {
    log2n = _log2n;
    length = 1 << log2n;
    increase = _increase;
    logDataNum = _logDataNum;
    dataNum = 1 << logDataNum;
    colNum = _colNum;
    setMaskNum();
}

MaskingGenerator::~MaskingGenerator() {
    if(mask != NULL) delete[] mask;
}

void MaskingGenerator::setMaskNum() {
    if (colNum == -1) {
        maskNum = log2n * (log2n + 1) / 2;
    } else {
        maskNum = (log2n - logDataNum) * (log2n - logDataNum + 1) / 2;
    }
    
    mask = new double*[maskNum];
    for(int i = 0; i < maskNum; i++) {
        mask[i] = new double[length]();
    }
}

void MaskingGenerator::printMask(double** mask, long maskNum) {
    for(int i = 0; i < maskNum; i++){
        std::cout << "mask[" << i << "] = [";
        for(int j = 0; j < length; j++) {
            std::cout << mask[i][j] << ", ";
        }
        std::cout << "]" << std::endl;        
    }
    std::cout << std::endl;
}

double** MaskingGenerator::getMasking() {
    generateMaskingRec(log2n, 0, 0);
    return mask;
}

double** MaskingGenerator::getMaskingOther() {
    generateMaskingRecOther(log2n, 0, 0);
    return mask;
}

double** MaskingGenerator::getBitonicMergeMasking() {
    for(int i = 0; i < log2n; i++) {
        generateBitonicMergeMasking(i);
    }
    return mask;
}

// double** MaskingGenerator::getTableMasking() {
//     long logNum =  log2n - logDataNum;
//     long maskNum = logNum * (logNum + 1) / 2;
//     setMaskNum(maskNum);
//     generateMaskingTableRec(logNum, 0, 0);
//     return mask;
// }

// double** MaskingGenerator::getTableMaskingBy(long colNum) {
//     long logNum =  log2n - logDataNum;
//     long maskNum = logNum * (logNum + 1) / 2;
//     setMaskNum(maskNum);
    
//     generateMaskingTableByRec(colNum, logNum, 0, 0);
//     return mask;
// }

void MaskingGenerator::generateBitonicMergeMasking(long num) {
    for(int j = 0; j < (1 << num); j++) {
        for(int k = 0; k < (1 << (log2n - 1 - num)); k++) {
            long loc = j * (1 << (log2n - num)) + k;
            if (!increase) loc = length - 1 - loc;
            mask[num][loc] = 1; 
        }
    }
}

void MaskingGenerator::generateMaskingComparison(long loc, long jump) {
    jump *= dataNum;
    long repeat = length / (jump * 2);
    
    for(int i = 0; i < repeat; i++) {
        for(int j = 0; j < jump; j++) {
            if (colNum == -1 || j % dataNum == colNum) {
                long tmp = i * jump * 2 + j;
                if(!increase) {
                    tmp = length - 1 - tmp;
                }
                mask[loc][tmp] = 1;
            }            
        }
    }
}

void MaskingGenerator::generateMaskingMerge(long loc, long num, long jump) {
    jump *= dataNum;
    long repeat = length / (jump * num);
    
    for(int i = 0; i < repeat; i++) {
        for(int j = 0; j < jump; j++) {
            if(colNum == -1 || j % dataNum == colNum) {
                for(int k = 0; k < num / 2 - 1; k++) {
                    long tmp = i * jump * num + (2 * k + 1) * jump + j;
                    if(!increase) {
                        tmp = length - 1 - tmp;
                    } 
                    mask[loc][tmp] = 1;    
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

long MaskingGenerator::generateMaskingRecOther(long logNum, long logJump, long loc) {
    if (logNum == 1) {
        generateMaskingComparisonOther(loc, 1 << logJump);
    } else {
        if (logJump == 0) {
            loc = generateMaskingRecOther(logNum - 1, logJump, loc);
        }
        loc = generateMaskingRecOther(logNum - 1, logJump + 1, loc);
        generateMaskingMergeOther(loc, 1 << logNum, 1 << logJump);
    }
    return loc + 1;
}

void MaskingGenerator::generateMaskingComparisonOther(long loc, long jump) {
    jump *= dataNum;
    long repeat = length / (jump * 2);
    
    for(int i = 0; i < repeat; i++) {
        for(int j = 0; j < jump; j++) {
            if (colNum == -1 || j % dataNum == colNum) {
                long tmp = i * jump * 2 + j + jump;
                if(!increase) {
                    tmp = length - 1 - tmp;
                }
                mask[loc][tmp] = 1;
            }            
        }
    }
}

void MaskingGenerator::generateMaskingMergeOther(long loc, long num, long jump) {
    jump *= dataNum;
    long repeat = length / (jump * num);
    
    for(int i = 0; i < repeat; i++) {
        for(int j = 0; j < jump; j++) {
            if(colNum == -1 || j % dataNum == colNum) {
                for(int k = 0; k < num / 2 - 1; k++) {
                    long tmp = i * jump * num + (2 * k + 1) * jump + j + jump;
                    if(!increase) {
                        tmp = length - 1 - tmp;
                    } 
                    mask[loc][tmp] = 1;    
                }            
            }
            
        }        
    }
}

// void MaskingGenerator::generateMaskingTableComp(long loc, long jump) {
//     generateMaskingComparison(loc, jump * (1 << logDataNum));
// }

// void MaskingGenerator::generateMaskingTableMerge(long loc, long num, long jump) {
//     generateMaskingMerge(loc, num, jump * (1 << logDataNum));
// }

// long MaskingGenerator::generateMaskingTableRec(long logNum, long logJump, long loc) {
//     if (logNum == 1) {
//         generateMaskingTableComp(loc, 1 << logJump);
//     } else {
//         if (logJump == 0) {
//             loc = generateMaskingTableRec(logNum - 1, logJump, loc);
//         }
//         loc = generateMaskingTableRec(logNum - 1, logJump + 1, loc);
//         generateMaskingTableMerge(loc, 1 << logNum, 1 << logJump);
//     }
//     return loc + 1;
// }

// long MaskingGenerator::generateMaskingTableByRec(long colNum, long logNum, long logJump, long loc) {
//     if (logNum == 1) {
//         generateMaskingTableCompBy(colNum, loc, 1 << logJump);
//     } else {
//         if (logJump == 0) {
//             loc = generateMaskingTableByRec(colNum, logNum - 1, logJump, loc);
//         }
//         loc = generateMaskingTableByRec(colNum, logNum - 1, logJump + 1, loc);
//         generateMaskingTableMergeBy(colNum, loc, 1 << logNum, 1 << logJump);
//     }
//     return loc + 1;
// }

// void MaskingGenerator::generateMaskingTableCompBy(long colNum, long loc, long jump) {
//     long dataNum = 1 << logDataNum;
//     jump *= dataNum;

//     long repeat = length / (jump * 2);
//     for(int i = 0; i < repeat; i++) {
//         for(int j = 0; j < jump; j++) {
//             if (j % dataNum == colNum) {
//                 long tmp = i * jump * 2 + j;
//                 if(!increase) {
//                     tmp = length - 1 - tmp;
//                 }
//                 mask[loc][tmp] = 1;
//             }
//         }
//     }
// }

// void MaskingGenerator::generateMaskingTableMergeBy(long colNum, long loc, long num, long jump) {
//     long dataNum = (1 << logDataNum);
//     jump *= dataNum;

//     long repeat = length / (jump * num);
//     for(int i = 0; i < repeat; i++) {
//         for(int j = 0; j < jump; j++) {
//             if (j % dataNum == colNum) {
//                 for(int k = 0; k < num / 2 - 1; k++) {
//                     long tmp = i * jump * num + (2 * k + 1) * jump + j;
//                     if(!increase) {
//                         tmp = length - 1 - tmp;
//                     } 
//                     mask[loc][tmp] = 1;    
//                 }  
//             }          
//         }        
//     }
// }