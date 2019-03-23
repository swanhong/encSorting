#include "test/TestPlain.h"
#include "test/TestEnc.h"
#include "test/TestBoot.h"
#include "test/TestSort.h"


int main() { 
    // ******************************
    // Parameters (long)
    // {logN, logQ, logp, logc, log2n, radix, logq, logT}
    // ******************************
    Parameter sortingTestParamSmall = {11, 1500, 30, 30, 6, 4, 35, 4};
    Parameter sortingTestParam1 = {12, 1000, 40, 40, 10, 32, 45, 5};
    Parameter sortingTestParamBig = {15, 1200, 40, 40, 14, 128, 45, 5};
    Parameter sortingTestParamBig2 = {16, 1200, 40, 40, 15, 32, 45, 4};
    Parameter sortingTestParam3 = {13, 2000, 30, 30, 12, 64, 35, 4};

    
    // ******************************
    // *** Test PlainSort
    // ******************************
    // TestPlain::showMasking(6);
    // TestPlain::plainSort(5);

    // ******************************
    // *** Test algorithms for encrypted data
    // ******************************
    // TestEnc::approxSqrt(sortingTestParam1, 10);
    // TestEnc::minMax(sortingTestParam1, 10);
    // TestEnc::compAndSwap(sortingTestParamSmall, 10);

    // ******************************
    // *** Test algorithms for encrypted data with Bootstrapping
    // ******************************
    // TestBoot::approxSqrt(sortingTestParam1, 12);
    // TestBoot::minMax(sortingTestParam1, 12);
    // TestBoot::compAndSwap(sortingTestParam1, 12);

    // ******************************
    // *** Test EncSorting
    // ******************************
    TestSort::sort(sortingTestParam1, 10);
    return 0;
}