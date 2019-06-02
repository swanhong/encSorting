#include "test/TestMask.h"
#include "test/TestPlain.h"
#include "test/TestAlgo.h"
#include "test/TestTableAlgo.h"
#include "test/TestSort.h"


int main() { 
    // ******************************
    // Parameters (long)
    // {logN, logQ, logp, logc, log2n, radix, logq, logT}
    // ******************************
    Parameter sortingTestParamSmall = {7, 1350, 40, 40, 4, 4, 50, 4};
    Parameter sortingTestParamSmall2 = {7, 1350, 40, 40, 6, 8, 50, 4};
    Parameter sortingTestParam1 = {12, 1350, 40, 40, 10, 32, 50, 5};
    Parameter sortingTestParamBig = {15, 650, 30, 30, 14, 128, 40, 5};
    Parameter sortingTestParamBig2 = {16, 1350, 40, 40, 15, 32, 50, 4};
    Parameter sortingTestParamBig3 = {17, 1350, 40, 40, 16, 256, 50, 4};

    // ******************************
    // *** Test Maskings
    // ******************************
    long logn = 5;
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
    // TestMask::showTableMasking(logn, logDataNum, true);
    // TestMask::showTableMaskingBy(logn, logDataNum, 0, true);
    // TestMask::showTableMaskingOther(logn + logDataNum, logDataNum, 0, true);
    // TestMask::showReverseMasking(logn, true);
    // TestMask::showReverseMaskingRight(logn, true);

    // TestMask::showReverseMasking(logn, false);
    // TestMask::showReverseMaskingRight(logn, false);

    // ******************************
    // *** Test PlainSort
    // ******************************
    
    // TestPlain::plainSort(5);
    // TestPlain::bitonicMerge(4, 4);
    // TestPlain::plainTableSort(4, 1, 0, true);
    // TestPlain::plainTableSort(4, 1, 0, false);

    // ******************************
    // *** Test algorithms for encrypted data
    // ******************************
    // ******************************
    // TestAlgo::approxSqrt(sortingTestParam1, 10);    
    // TestAlgo::minMax(sortingTestParamSmall, 15);
    // TestAlgo::EncSwap(sortingTestParamSmall, 13);
    // TestAlgo::reverse(sortingTestParamSmall);
    // TestAlgo::halfCleaner(sortingTestParamSmall, 15);
    // TestAlgo::approxInverse(sortingTestParamSmall, 5);
    // TestAlgo::comparison(sortingTestParamSmall, 10, 5);
    // TestAlgo::encSwapTable(sortingTestParamSmall, 2, 0, 5, 5);

    // ******************************
    // *** Check Parameters
    // ******************************
    // TestAlgo::bootstrapping(sortingTestParamSmall);
    // TestEnc::compAndSwap(sortingTestParam1, 10);

    // ******************************
    // *** Test EncSorting
    // ******************************
    // TestSort::sort(sortingTestParamSmall, 10);
    TestSort::merge(sortingTestParamSmall, 13, 2);
    // TestSort::sortAndMerge(sortingTestParamSmall, 15, 4);

    // TestSort::tableSort(sortingTestParamSmall, 1, 0, 12, 5); 
    // TestSort::tableMerge(sortingTestParamSmall, 2, 2, 0, 5, 5);    
    // TestSort::bitonicSort(sortingTestParamSmall, 13);
    
    return 0;
}