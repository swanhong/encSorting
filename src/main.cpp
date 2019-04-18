#include "test/TestPlain.h"
#include "test/TestEnc.h"
#include "test/TestBoot.h"
#include "test/TestSort.h"

int main() { 
    // ******************************
    // Parameters (long)
    // {logN, logQ, logp, logc, log2n, radix, logq, logT}
    // ******************************
    Parameter sortingTestParamSmall = {6, 3000, 40, 40, 2, 2, 45, 4};
    Parameter sortingTestParam1 = {12, 1000, 40, 40, 10, 32, 45, 5};
    Parameter sortingTestParamBig = {15, 1200, 40, 40, 14, 128, 45, 5};
    Parameter sortingTestParamBig2 = {16, 1200, 40, 40, 15, 32, 45, 4};
    Parameter sortingTestParam3 = {13, 2000, 30, 30, 12, 64, 35, 4};

    
    // ******************************
    // *** Test PlainSort
    // ******************************
    // TestPlain::showMasking(5, true);
    // TestPlain::showMasking(5, false);
    // TestPlain::showBitonicMergeMasking(6);
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
    // TestBoot::minMax(sortingTestParam1, 12);
    // TestBoot::compAndSwap(sortingTestParamSmall, 15);
    // TestBoot::reverse(sortingTestParamSmall);

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
    return 0;
}