#include "test/TestMask.h"
#include "test/TestPlain.h"
#include "test/TestEnc.h"
#include "test/TestBoot.h"
#include "test/TestSort.h"


int main() { 
    // ******************************
    // Parameters (long)
    // {logN, logQ, logp, logc, log2n, radix, logq, logT}
    // ******************************
    Parameter sortingTestParamSmall = {7, 1350, 40, 40, 4, 4, 50, 4};
    Parameter sortingTestParamSmall2 = {7, 1350, 50, 50, 6, 8, 65, 4};
    Parameter paramLee = {12, 1350, 50, 50, 10, 32, 60, 5};
    Parameter sortingTestParam1 = {10, 1350, 45, 45, 6, 8, 55, 5};
    Parameter sortingTestParamBig = {15, 1350, 40, 40, 14, 128, 50, 5};
    Parameter sortingTestParamBig2 = {16, 1350, 40, 40, 15, 32, 50, 4};
    Parameter sortingTestParamBig3 = {17, 1350, 40, 40, 16, 256, 50, 4};


    // ******************************
    // *** Test Maskings
    // ******************************
    long logn = 4;
    long logDataNum = 2;
    long colNum = 1;
    // TestMask::showMasking(logn, true);
    // TestMask::showMasking(logn, false);
    // TestMask::showMaskingOther(logn, true);
    // TestMask::showBitonicMergeMasking(logn, true);
    // TestMask::showBitonicMergeMasking(logn, false);
    // TestMask::showColNumMasking(logn, logDataNum, colNum, true);
    // TestMask::showTableMergeMasking(logn, logDataNum, colNum, true);
    // TestMask::showTableMergeMaskingOther(logn, logDataNum, colNum, true);
    // TestMask::showTableMergeMasking(logn, logDataNum, colNum, false);
    // TestMask::showTableMergeMaskingOther(logn, logDataNum, colNum, false);
    // TestMask::showTableMasking(logn + logDataNum, logDataNum, true);
    // TestMask::showTableMaskingBy(logn + logDataNum, logDataNum, 0, true);
    // TestMask::showTableMaskingOther(logn + logDataNum, logDataNum, 0, true);

    // ******************************
    // *** Test PlainSort
    // ******************************
    
    // TestPlain::plainSort(5);
    // TestPlain::bitonicMerge(4, 4);
    // TestPlain::plainTableSort(4, 1, 0, true);
    // TestPlain::plainTableSort(4, 1, 0, false);

    // ******************************
    // *** Test algorithms for encrypted data with Bootstrapping
    // ******************************
    // ******************************
    // TestBoot::approxSqrt(sortingTestParam1, 15);
    // TestBoot::approxInverse(sortingTestParamSmall, 5);
    // TestBoot::approxComp(sortingTestParamSmall, 5, 5);
    // TestBoot::minMax(sortingTestParamSmall, 15);
    // TestBoot::compAndSwap(sortingTestParamSmall, 13);
    // TestBoot::reverse(sortingTestParamSmall);
    // TestBoot::compAndSwapTable(sortingTestParamSmall, 2, 0, 5, 5);

    // ******************************
    // *** Check Parameters
    // ******************************
    // TestBoot::bootstrapping(sortingTestParam1);
    // TestEnc::compAndSwap(sortingTestParam1, 10);

    // ******************************
    // *** Test EncSorting
    // ******************************
    TestSort::sort(sortingTestParam1, 13);
    // TestSort::merge(sortingTestParamSmall, 13, 2);
    // TestSort::sortAndMerge(sortingTestParamSmall, 15, 4);

    // TestSort::tableSort(sortingTestParam1, 0, 0, 4, 4, true);
    // TestSort::tableMerge(sortingTestParamSmall, 2, 2, 0, 5, 5);    

    return 0;
}