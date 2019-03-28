#ifndef TESTSORT_H_
#define TESTSORT_H_

#include "../enc/EncSorting.h"
#include "../plain/PlainSorting.h"

class TestSort {
public:
    static void sort(Parameter param, long iter, bool=true);
    static void bitonicMerge(Parameter param, long iter);
    static void testMerge(Parameter param, long iter, long logNum);
};



#endif // !TESTSORT_H_