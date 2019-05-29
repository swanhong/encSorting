#ifndef TestAlgo_H_
#define TestAlgo_H_

#include "../Parameter.h"
#include "../PrintUtils.h"
#include "../enc/EncAlgo.h"
#include "../plain/PlainSorting.h"
#include "../MaskingGenerator.h"
#include "stdlib.h"


class TestAlgo {
public:
    static void bootstrapping(Parameter parameter);
	static void approxSqrt(Parameter parameter, long iter);
	static void minMax(Parameter parameter, long iter);
	static void EncSwap(Parameter param, long iter);
	static void reverse(Parameter param);
	static void halfCleaner(Parameter param, long iter);

	//************
	
	static void approxInverse(Parameter param, long iter);
	static void comparison(Parameter param, long invIter, long compIter);
	static void encSwapTable(Parameter param, long logDataNum, long colNum, long invIter, long compIter);

};

#endif // !TestAlgo_H_