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
    
    Parameter sortingTestParam = {12, 850, 30, 30, 4, 4, 35, 4};

    // TestBootstrapping::testSqrtWithBoot(sortingTestParam, 20);
    // TestBootstrapping::testMaxMinWithBoot(sortingTestParam, 20);
    // TestBootstrapping::testMaxMinWithBootAndDecrypt(bootstrapping_test_param4, 10);
    // TestBootstrapping::testEncCompAndSwapWithBoot(sortingTestParam, 20);
    

    // TestBootstrapping::bootstrapping_test(sortingTestParam);
    TestBootstrapping::bootstrapping_test_with_mult(sortingTestParam);

    // TestBootstrapping::testEncSort(sortingTestParam, 5);
    // TestBootstrapping::testEncSortWithDecrypt(sortingTestParam, 3);



    return 0;
}