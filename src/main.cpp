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
    
    Parameter sortingTestParamSmall = {8, 850, 40, 40, 6, 8, 45, 4};
    Parameter sortingTestParam1 = {12, 2000, 30, 30, 4, 16, 35, 4};
    Parameter sortingTestParamBig = {15, 1200, 30, 30, 14, 128, 35, 4};
    Parameter sortingTestParam3 = {13, 2000, 30, 30, 12, 64, 35, 4};

    // TestBootstrapping::testSqrtWithBoot(sortingTestParam, 20);
    // TestBootstrapping::testMaxMinWithBoot(sortingTestParam, 20);
    
    // TestBootstrapping::testEncCompAndSwapWithBoot(sortingTestParam, 20);
    

    // TestBootstrapping::bootstrapping_test(sortingTestParamBig);
    // TestBootstrapping::bootstrapping_test_with_mult(sortingTestParam1, 30);

    // TestBootstrapping::testEncSort(sortingTestParam1, 20);

    // TestBootstrapping::testSqrtWithBootAndDecrypt(sortingTestParam1, 30);
    // TestBootstrapping::testMaxMinWithBootAndDecrypt(sortingTestParam1, 10);
    TestBootstrapping::testEncSortWithDecrypt(sortingTestParam1, 8);
    // TestBootstrapping::testEncCompAndSwapWithBootAndDecrypt(sortingTestParamBig, 9);




    return 0;
}