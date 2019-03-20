#ifndef TESTBOOTSTRAPPING_H_
#define TESTBOOTSTRAPPING_H_

#include "Parameter.h"
#include "EncContext.h"
#include "BootContext.h"
#include "stdlib.h"
#include "SortingAlgorithm.h"
#include "TestAlgorithm.h"

class TestBootstrapping {
public:
    static void bootstrapping_test(Parameter parameter);

	static void bootstrapping_test_with_mult(Parameter parameter, int iter);
	
	static void testSqrtWithBoot(Parameter parameter, long iter);
	
	static void testMaxMinWithBoot(Parameter parameter, long iter);

	static void testEncCompAndSwapWithBoot(Parameter parameter, long iter);

	static void testEncSort(Parameter parameter, long iter);

	static void testEncSortWithDecrypt(Parameter parameter, long iter);

	static void testMaxMinWithBootAndDecrypt(Parameter parameter, long iter);

	static void testSqrtWithBootAndDecrypt(Parameter parameter, long iter);

	static void testEncCompAndSwapWithBootAndDecrypt(Parameter parameter, long iter);
};

#endif // !TESTBOOTSTRAPPING_H_