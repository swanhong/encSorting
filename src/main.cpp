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
    Parameter sortingTestParamSmall = {6, 1000, 40, 40, 4, 4, 45, 4};
    Parameter sortingTestParam1 = {12, 1000, 40, 40, 10, 32, 45, 5};
    Parameter sortingTestParamBig = {15, 1200, 40, 40, 14, 128, 45, 5};
    Parameter sortingTestParamBig2 = {16, 1200, 40, 40, 15, 32, 45, 4};
    Parameter sortingTestParamBig3 = {17, 1200, 40, 40, 16, 256, 45, 4};
    Parameter sortingTestParam3 = {13, 2000, 30, 30, 12, 64, 35, 4};


    // ******************************
    // *** Test Maskings
    // ******************************
    // long logn = 3;
    // long logDataNum = 2;
    // TestMask::showMasking(logn, true);
    // TestMask::showMasking(logn, false);
    // TestMask::showMaskingOther(logn, true);
    // TestMask::showBitonicMergeMasking(logn);
    // TestMask::showTableMasking(logn + logDataNum, logDataNum, true);
    // TestMask::showTableMaskingBy(logn + logDataNum, logDataNum, 0, true);
    // TestMask::showTableMaskingOther(logn + logDataNum, logDataNum, 0, true);

    // ******************************
    // *** Test PlainSort
    // ******************************
    
    // TestPlain::plainSort(5);
    // TestPlain::bitonicMerge(4, 4);

    // ******************************
    // *** Test algorithms for encrypted data
    // ******************************
    // TestEnc::approxSqrt(sortingTestParam1, 10);
    // TestEnc::minMax(sortingTestParam1, 10);
    // TestEnc::compAndSwap(sortingTestParamSmall, 10);

    // ******************************
    // *** Test algorithms for encrypted data with Bootstrapping
    // ******************************
    // ******************************
    // TestBoot::approxSqrt(sortingTestParamSmall, 15);
    // TestBoot::approxInverse(sortingTestParamSmall, 3);
    // TestBoot::approxComp(sortingTestParamSmall, 3, 6);
    // TestBoot::minMax(sortingTestParam1, 12);
    // TestBoot::compAndSwap(sortingTestParamSmall, 15);
    // TestBoot::reverse(sortingTestParamSmall);
    TestBoot::compAndSwapTable(sortingTestParam1, 2, 0, 5, 6);

    // ******************************
    // *** Check Parameters
    // ******************************
    // TestBoot::bootstrapping(sortingTestParamBig3);
    // TestEnc::compAndSwap(sortingTestParam1, 10);

    // ******************************
    // *** Test EncSorting
    // ******************************
    // TestSort::sort(sortingTestParamSmall, 15);
    // TestSort::merge(sortingTestParamSmall, 15, 4);
    // TestSort::sortAndMerge(sortingTestParamSmall, 15, 4);
    // TestSort::tableSort(sortingTestParamSmall, 2, 0, 5, 6);

    return 0;
}