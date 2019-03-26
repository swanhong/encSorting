#ifndef TestBoot_H_
#define TestBoot_H_

#include "../Parameter.h"
#include "../PrintUtils.h"
#include "../enc/BootAlgorithm.h"
#include "../MaskingGenerator.h"
#include "stdlib.h"


class TestBoot {
public:
    static void bootstrapping(Parameter parameter);
	
	static void approxSqrt(Parameter parameter, long iter);
	
	static void minMax(Parameter parameter, long iter);

	static void compAndSwap(Parameter param, long iter);
};

#endif // !TestBoot_H_