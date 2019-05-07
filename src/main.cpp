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
    Parameter sortingTestParamSmall = {7, 1350, 50, 50, 4, 4, 60, 4};
    Parameter sortingTestParamSmall2 = {7, 1350, 40, 40, 6, 8, 50, 4};
    Parameter paramLee = {12, 1350, 50, 50, 10, 32, 60, 5};
    Parameter sortingTestParam1 = {10, 1350, 40, 40, 6, 16, 50, 5};
    Parameter sortingTestParamBig = {15, 1350, 40, 40, 14, 128, 50, 5};
    Parameter sortingTestParamBig2 = {16, 1350, 40, 40, 15, 32, 50, 4};
    Parameter sortingTestParamBig3 = {17, 1350, 40, 40, 16, 256, 50, 4};


    // ******************************
    // *** Test Maskings
    // ******************************
    // long logn = 5;
    // long logDataNum = 2;
    // TestMask::showMasking(logn, true);
    // TestMask::showMasking(logn, false);
    // TestMask::showMaskingOther(logn, true);
    // TestMask::showBitonicMergeMasking(logn, true);
    // TestMask::showBitonicMergeMasking(logn, false);
    // TestMask::showTableMergeMasking(logn, logDataNum, 0, true);
    // TestMask::showTableMergeMasking(logn, logDataNum, 0, false);
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
    // TestBoot::approxSqrt(sortingTestParamSmall2, 15);
    // TestBoot::approxInverse(sortingTestParamSmall, 5);
    // TestBoot::approxComp(sortingTestParamSmall, 5, 5);
    // TestBoot::minMax(sortingTestParamSmall, 15);
    // TestBoot::compAndSwap(sortingTestParamSmall, 13);
    // TestBoot::reverse(sortingTestParamSmall);
    // TestBoot::compAndSwapTable(sortingTestParamSmall, 2, 0, 5, 5);

    // ******************************
    // *** Check Parameters
    // ******************************
    // TestBoot::bootstrapping(param);
    // TestEnc::compAndSwap(sortingTestParam1, 10);

    // ******************************
    // *** Test EncSorting
    // ******************************
    // TestSort::sort(sortingTestParamSmall, 13);
    // TestSort::merge(sortingTestParamSmall, 13, 2);
    // TestSort::sortAndMerge(sortingTestParamSmall, 15, 4);
    TestSort::tableSort(sortingTestParamSmall, 2, 2, 5, 5, false);

    return 0;
}