#include "TestAlgorithm.h"
#include "TestBootstrapping.h"

int main() {
    // TestAlgorithm::testMult(1000, 30, 15);
    
    // TestAlgorithm::testMasking(3);
    
    // TestAlgorithm::testPlainSort(5);
    
    // TestAlgorithm::testSqrt(bootstrapping_test_param4);

    // TestAlgorithm::testMaxMin(bootstrapping_test_param4);
    
    // TestAlgorithm::testEncCompAndSwap(bootstrapping_test_param4);

    // TestAlgorithm::testEncSorting(bootstrapping_test_param4);

    /* 
        Parameters (long)
        {logN, logQ, logp, logc, log2n, radix, logq, logT}
    */

    Parameter sortingTestParamSmall = {11, 1500, 30, 30, 6, 4, 35, 4};
    Parameter sortingTestParam1 = {12, 1000, 40, 40, 10, 32, 45, 5};
    Parameter sortingTestParamBig = {15, 1200, 40, 40, 14, 128, 45, 5};
    Parameter sortingTestParamBig2 = {16, 1200, 40, 40, 15, 32, 45, 4};

    Parameter sortingTestParam3 = {13, 2000, 30, 30, 12, 64, 35, 4};

    // TestBootstrapping::testSqrtWithBoot(sortingTestParam, 20);
    // TestBootstrapping::testMaxMinWithBoot(sortingTestParam, 20);
    
    // TestBootstrapping::testEncCompAndSwapWithBoot(sortingTestParam, 20);  

    // TestBootstrapping::bootstrapping_test(sortingTestParamBig2);
    // TestAlgorithm::testEncCompAndSwap(sortingTestParamBig2, 10);
    // TestBootstrapping::testMaxMinWithBootAndDecrypt(sortingTestParamSmall, 20);
    // TestBootstrapping::bootstrapping_test_with_mult(sortingTestParam1, 30);

    // TestBootstrapping::testEncSort(sortingTestParamSmall, 10);

    // TestBootstrapping::testSqrtWithBootAndDecrypt(sortingTestParamSmall, 10);
    
    TestBootstrapping::testEncSortWithDecrypt(sortingTestParamSmall, 10);
    
    // TestBootstrapping::testEncCompAndSwapWithBootAndDecrypt(sortingTestParamBig, 9);




    return 0;
}