#include "TestAlgorithm.h"
#include "TestBootstrapping.h"

// #include "../bootsrc/new_bootstrapping.h"
// #include "Parameter.h"

int main() {
    // TestAlgorithm::testMult(1000, 30, 15);
    
    // TestAlgorithm::testMasking(3);
    
    // TestAlgorithm::testPlainSort(5);
    
    // TestAlgorithm::testSqrt(bootstrapping_test_param4);

    // TestAlgorithm::testMaxMin(bootstrapping_test_param4);
    
    // TestAlgorithm::testEncCompAndSwap(bootstrapping_test_param4);

    // TestAlgorithm::testEncSorting(bootstrapping_test_param4);
    
    Parameter sortingTestParamSmall = {11, 1500, 30, 30, 6, 4, 35, 4};
    Parameter sortingTestParam1 = {12, 2000, 30, 30, 4, 16, 35, 4};
    Parameter sortingTestParamBig = {15, 1200, 40, 40, 14, 128, 45, 4};
    Parameter sortingTestParamBig2 = {16, 1200, 30, 30, 15, 32, 35, 4};

    Parameter sortingTestParam3 = {13, 2000, 30, 30, 12, 64, 35, 4};

    // TestBootstrapping::testSqrtWithBoot(sortingTestParam, 20);
    // TestBootstrapping::testMaxMinWithBoot(sortingTestParam, 20);
    
    // TestBootstrapping::testEncCompAndSwapWithBoot(sortingTestParam, 20);
    

    // TestBootstrapping::bootstrapping_test(sortingTestParamSmall);
    // TestBootstrapping::bootstrapping_test_with_mult(sortingTestParam1, 30);

    // TestBootstrapping::testEncSort(sortingTestParam1, 20);

    // TestBootstrapping::testSqrtWithBootAndDecrypt(sortingTestParam1, 30);
    // TestBootstrapping::testMaxMinWithBootAndDecrypt(sortingTestParam1, 10);
    TestBootstrapping::testEncSortWithDecrypt(sortingTestParamSmall, 20);
    // TestAlgorithm::testEncCompAndSwap(sortingTestParamSmall, 12);
    // TestBootstrapping::testEncCompAndSwapWithBootAndDecrypt(sortingTestParamBig, 9);




    return 0;
}