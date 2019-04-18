#ifndef TESTSORT_H_
#define TESTSORT_H_

#include "../enc/EncSorting.h"
#include "../plain/PlainSorting.h"

class TestSort {
public:
    static void sort(Parameter param, long iter, bool=true);
    static void merge(Parameter param, long iter, long logNum);
    static void sortAndMerge(Parameter param, long iter, long logNum);
};



#endif // !TESTSORT_H_