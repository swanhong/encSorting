#ifndef TestTableAlgo_H_
#define TestTableAlgo_H_

#include "../Parameter.h"
#include "../PrintUtils.h"
#include "../enc/EncTableAlgo.h"
#include "../MaskingGenerator.h"
#include "stdlib.h"


class TestTableAlgo {
public:
	static void approxInverse(Parameter parameter, long iter);

	static void approxComp(Parameter parameter, long invIter, long compIter);
	
	static void compAndSwapTable(Parameter param, long logDataNum, long colNum, long invIter, long compIter);

	static void runTestAlgo(Parameter param, long* iter);
};

#endif // !TestAlgo_H_