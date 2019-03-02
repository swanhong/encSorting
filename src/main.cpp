// #include "iostream"
#include "TestAlgorithm.h"
// #include "improved_bootstrapping/new_bootstrapping.h"
// #include "TestLogicHeaan.h"


// #include "LogicHeaan.h"

using namespace std;

Parameter bootstrapping_test_param1 = {16, 850, 30, 30, 6, 4, 35, 4};
Parameter bootstrapping_test_param2 = {16, 850, 30, 30, 9, 8, 35, 4};
Parameter bootstrapping_test_param3 = {16, 850, 30, 30, 12, 16, 35, 4};
Parameter bootstrapping_test_param4 = {12, 800, 30, 30, 4, 16, 35, 4};

int main() {
    
    // TestAlgorithm::testMult(1000, 30, 15);
    
    // TestAlgorithm::testMasking(3);
    
    TestAlgorithm::testPlainSort(5);
    
    // TestAlgorithm::testSqrt(bootstrapping_test_param4);

    // TestAlgorithm::testMaxMin(bootstrapping_test_param4);
    
    // TestAlgorithm::testEncCompAndSwap(bootstrapping_test_param4);

    // TestAlgorithm::testEncSorting(bootstrapping_test_param4);


    return 0;
}



