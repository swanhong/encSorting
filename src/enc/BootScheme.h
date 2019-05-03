#ifndef BOOTSCHEME_H_
#define BOOTSCHEME_H_

#include "../../HEAAN/src/HEAAN.h"
#include "../bootsrc/new_bootstrapping.h"
#include "../Parameter.h"
#include "../PrintUtils.h"

#include <string>

class BootScheme : public Scheme {
private:
    long NUM_OF_MULT = 0;
    long NUM_OF_BOOT = 0;
    long NUM_OF_CURRENT_MULT = 0;
    long NUM_OF_CURRENT_BOOT = 0;

public:
    BootScheme(SecretKey& secretKey, Ring& ring) : Scheme(secretKey, ring) {}
    ~BootScheme() {}

    void countMult();
    void countBoot();
    void resetCount();

    void checkAndBoot(Ciphertext& cipher, bool condition, BootHelper& bootHelper, Parameter param);

    Ciphertext multWithBoot(Ciphertext& cipher1, Ciphertext& cipher2, BootHelper& bootHelper, Parameter param);
    void multAndEqualWithBoot(Ciphertext& cipher1, Ciphertext& cipher2, BootHelper& bootHelper, Parameter param);
    
    void modDownToAndEqualModified(Ciphertext& cipher1, Ciphertext& cipher2, BootHelper& bootHelper, Parameter param);
    
    Ciphertext squareWithBoot(Ciphertext& cipher, BootHelper& bootHelper, Parameter param);
    void squareAndEuqalWithBoot(Ciphertext& cipher, BootHelper& bootHelper, Parameter param);

    Ciphertext multByPolyWithBoot(Ciphertext& cipher, ZZ* poly, BootHelper& bootHelper, Parameter param);
    void multByPolyAndEqualWithBoot(Ciphertext& cipher, ZZ* poly, BootHelper& bootHelper, Parameter param);

    void nomalizeAndEuqal(Ciphertext& cipher);

    void decryptAndPrint(std::string str, SecretKey& secretKey, Ciphertext& cipher);

    void showTotalCount();
    void showCurrentCount();
};

#endif // !BOOTSCHEME_H_
