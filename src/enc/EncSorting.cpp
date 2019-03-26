#include "EncSorting.h"


void EncSorting::runSorting(Ciphertext& cipher, Parameter param, long iter, BootScheme& scheme, Ring& ring, BootHelper& bootHelper) {
        MaskingGenerator mg(param.log2n);
        double** mask = mg.getMasking();
        BootAlgo bootAlgo;
        sortingRecursion(cipher, param.log2n, 0, 0, param, iter, mask, bootAlgo, scheme, ring, bootHelper);
        scheme.showTotalCount();
}

long EncSorting::sortingRecursion(Ciphertext& cipher, long logNum, long logJump, long loc, Parameter param, long iter, double** mask, BootAlgo& bootAlgo, BootScheme& scheme, Ring& ring, BootHelper& bootHelper) {
    TimeUtils timeutils;
    PrintUtils::nprint("run encBatcherSort with loc = " + to_string(loc), WANT_TO_PRINT);
    if (logNum == 1) {
        timeutils.start(to_string(loc)+"th CompAndSwap");
        bootAlgo.compAndSwap(cipher, mask[loc], 1<< logJump, iter, param, scheme, ring, bootHelper);
        scheme.showCurrentCount();
        scheme.resetCount();
        timeutils.stop(to_string(loc)+"th CompAndSwap");
    } else {
        if (logJump == 0) {
            loc = sortingRecursion(cipher, logNum - 1, logJump, loc, param, iter, mask, bootAlgo, scheme, ring, bootHelper);
        }
        loc = sortingRecursion(cipher, logNum - 1, logJump + 1, loc, param, iter, mask, bootAlgo, scheme, ring, bootHelper);

        timeutils.start(to_string(loc)+"th CompAndSwap");
        bootAlgo.compAndSwap(cipher, mask[loc], 1<< logJump, iter, param, scheme, ring, bootHelper);
        scheme.showCurrentCount();
        scheme.resetCount();
        timeutils.stop(to_string(loc)+"th CompAndSwap");
    }
    return loc + 1;
}