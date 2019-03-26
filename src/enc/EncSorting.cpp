#include "EncSorting.h"


void EncSorting::runSorting(Ciphertext& sortedCipher, const Ciphertext& inCipher, Parameter param, long iter, BootScheme& scheme, Ring& ring, BootHelper& bootHelper) {
        MaskingGenerator mg(param.log2n);
        double** mask = mg.getMasking();
        BootAlgo bootAlgo;
        sortingRecursion(sortedCipher, inCipher, param.log2n, 0, 0, param, iter, mask, bootAlgo, scheme, ring, bootHelper);
        scheme.showTotalCount();
}

long EncSorting::sortingRecursion(Ciphertext& sortedCipher, const Ciphertext& inCipher, long logNum, long logJump, long loc, Parameter param, long iter, double** mask, BootAlgo& bootAlgo, BootScheme& scheme, Ring& ring, BootHelper& bootHelper) {
    TimeUtils timeutils;
    Ciphertext a = inCipher;
    Ciphertext dummy;
    PrintUtils::nprint("run encBatcherSort with loc = " + to_string(loc), WANT_TO_PRINT);
    if (logNum == 1) {
        timeutils.start(to_string(loc)+"th CompAndSwap");
        bootAlgo.compAndSwap(dummy, a, mask[loc], 1<< logJump, iter, param, scheme, ring, bootHelper);
        scheme.showCurrentCount();
        scheme.resetCount();
        timeutils.stop(to_string(loc)+"th CompAndSwap");
        a = dummy;
    } else {
        if (logJump == 0) {
            loc = sortingRecursion(dummy, a, logNum - 1, logJump, loc, param, iter, mask, bootAlgo, scheme, ring, bootHelper);
            a = dummy;
        }
        loc = sortingRecursion(dummy, a, logNum - 1, logJump + 1, loc, param, iter, mask, bootAlgo, scheme, ring, bootHelper);
        a = dummy;

        timeutils.start(to_string(loc)+"th CompAndSwap");
        bootAlgo.compAndSwap(dummy, a, mask[loc], 1<< logJump, iter, param, scheme, ring, bootHelper);
        scheme.showCurrentCount();
        scheme.resetCount();
        timeutils.stop(to_string(loc)+"th CompAndSwap");
        a = dummy;
    }
    sortedCipher = a;
    return loc + 1;
}