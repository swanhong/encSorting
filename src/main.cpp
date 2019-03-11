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
    
    Parameter sortingTestParamSmall = {8, 850, 40, 30, 6, 8, 45, 4};
    Parameter sortingTestParam1 = {12, 800, 30, 30, 4, 16, 35, 4};
    Parameter sortingTestParamBig = {16, 1500, 30, 30, 15, 32, 35, 4};
    Parameter sortingTestParam3 = {13, 800, 30, 30, 12, 64, 35, 4};

    // TestBootstrapping::testSqrtWithBoot(sortingTestParam, 20);
    // TestBootstrapping::testMaxMinWithBoot(sortingTestParam, 20);
    
    // TestBootstrapping::testEncCompAndSwapWithBoot(sortingTestParam, 20);
    

    // TestBootstrapping::bootstrapping_test(sortingTestParamSmall);
    // TestBootstrapping::bootstrapping_test_with_mult(sortingTestParamSmall, 20);

    // TestBootstrapping::testEncSort(sortingTestParam1, 20);

    // TestBootstrapping::testSqrtWithBootAndDecrypt(sortingTestParamSmall, 30);
    // TestBootstrapping::testMaxMinWithBootAndDecrypt(sortingTestParamSmall, 30);
    // TestBootstrapping::testEncSortWithDecrypt(sortingTestParamSmall, 1);
    TestBootstrapping::testEncCompAndSwapWithBootAndDecrypt(sortingTestParamSmall, 30);




    return 0;
}