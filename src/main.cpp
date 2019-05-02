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
    Parameter sortingTestParamSmall = {6, 1350, 40, 40, 4, 4, 50, 4};
    Parameter sortingTestParam1 = {12, 1350, 40, 40, 8, 16, 50, 5};
    Parameter sortingTestParamBig = {15, 1350, 40, 40, 14, 128, 50, 5};
    Parameter sortingTestParamBig2 = {16, 1350, 40, 40, 15, 32, 50, 4};
    Parameter sortingTestParamBig3 = {17, 1350, 40, 40, 16, 256, 50, 4};


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
    // TestPlain::plainTableSort(4, 1);

    // ******************************
    // *** Test algorithms for encrypted data with Bootstrapping
    // ******************************
    // ******************************
    // TestBoot::approxSqrt(sortingTestParam1, 25);
    // TestBoot::approxInverse(sortingTestParamSmall, 5);
    // TestBoot::approxComp(sortingTestParamSmall, 5, 5);
    // TestBoot::minMax(sortingTestParamSmall, 20);
    TestBoot::compAndSwap(sortingTestParam1, 15);
    // TestBoot::reverse(sortingTestParamSmall);
    // TestBoot::compAndSwapTable(sortingTestParamSmall, 2, 0, 4, 5);

    // ******************************
    // *** Check Parameters
    // ******************************
    // TestBoot::bootstrapping(sortingTestParamSmall);
    // TestEnc::compAndSwap(sortingTestParam1, 10);

    // ******************************
    // *** Test EncSorting
    // ******************************
    // TestSort::sort(sortingTestParamSmall, 15);
    // TestSort::merge(sortingTestParamSmall, 15, 4);
    // TestSort::sortAndMerge(sortingTestParamSmall, 15, 4);
    // TestSort::tableSort(sortingTestParam1, 4, 0, 5, 5);

    return 0;
}