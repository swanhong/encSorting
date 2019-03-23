#ifndef BOOTSCHEME_H_
#define BOOTSCHEME_H_

#include "../../HEAAN/src/HEAAN.h"
#include "../bootsrc/new_bootstrapping.h"
#include "../Parameter.h"
#include "../PrintUtils.h"

class BootScheme : public Scheme {
public:
    BootScheme(SecretKey& secretKey, Ring& ring) : Scheme(secretKey, ring) {}
    ~BootScheme() {}

    void checkAndBoot(Ciphertext& cipher, bool condition, BootHelper& bootHelper, Parameter param);
    void multAndEqualWithBoot(Ciphertext& cipher1, Ciphertext& cipher2, BootHelper& bootHelper, Parameter param);
    void modDownToAndEqualModified(Ciphertext& cipher1, Ciphertext& cipher2, BootHelper& bootHelper, Parameter param);
    void squareAndEuqalWithBoot(Ciphertext& cipher, BootHelper& bootHelper, Parameter param);
    void multByPolyAndEqualWithBoot(Ciphertext& cipher, ZZ* poly, BootHelper& bootHelper, Parameter param);
};

#endif // !BOOTSCHEME_H_
