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

    // TestBootstrapping::testSqrtWithBoot(bootstrapping_test_param4, 20);
    // TestBootstrapping::testMaxMinWithBoot(bootstrapping_test_param4, 10);
    // TestBootstrapping::testEncCompAndSwapWithBoot(bootstrapping_test_param4, 10);
    
    
    Parameter sortingTestParam = {12, 1600, 30, 30, 9, 8, 35, 4};

    // TestBootstrapping::bootstrapping_test(sortingTestParam);

    // TestBootstrapping::testEncSort(sortingTestParam, 9);
    TestBootstrapping::testEncSortWithDecrypt(sortingTestParam, 12);



    return 0;
}



