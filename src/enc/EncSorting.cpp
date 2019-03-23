#include "EncSorting.h"


void EncSorting::runSorting(Ciphertext& sortedCipher, const Ciphertext& inCipher, Parameter param, long iter, BootScheme& scheme, Ring& ring, BootHelper& bootHelper) {
        MaskingGenerator mg(param.log2n);
        double** mask = mg.getMasking();
        BootAlgo bootAlgo;
        sortingRecursion(sortedCipher, inCipher, param.log2n, 0, 0, param, iter, mask, bootAlgo, scheme, ring, bootHelper);
}

long EncSorting::sortingRecursion(Ciphertext& sortedCipher, const Ciphertext& inCipher, long logNum, long logJump, long loc, Parameter param, long iter, double** mask, BootAlgo& bootAlgo, BootScheme& scheme, Ring& ring, BootHelper& bootHelper) {
    Ciphertext a = inCipher;
    Ciphertext dummy;
    PrintUtils::nprint("run encBatcherSort with loc = " + to_string(loc), WANT_TO_PRINT);
    if (logNum == 1) {
        bootAlgo.compAndSwap(dummy, a, mask[loc], 1<< logJump, iter, param, scheme, ring, bootHelper);
        a = dummy;
    } else {
        if (logJump == 0) {
            loc = sortingRecursion(dummy, a, logNum - 1, logJump, loc, param, iter, mask, bootAlgo, scheme, ring, bootHelper);
            a = dummy;
        }
        loc = sortingRecursion(dummy, a, logNum - 1, logJump + 1, loc, param, iter, mask, bootAlgo, scheme, ring, bootHelper);
        a = dummy;

        bootAlgo.compAndSwap(dummy, a, mask[loc], 1<< logJump, iter, param, scheme, ring, bootHelper);
        a = dummy;
    }
    sortedCipher = a;
    return loc + 1;
}