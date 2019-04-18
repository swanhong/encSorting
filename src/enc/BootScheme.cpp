#include "BootScheme.h"

void BootScheme::countMult() {
    NUM_OF_MULT += 1;
    NUM_OF_CURRENT_MULT += 1;
}

void BootScheme::countBoot() {
    NUM_OF_BOOT += 1;
    NUM_OF_CURRENT_BOOT += 1;
}

void BootScheme::resetCount() {
    NUM_OF_CURRENT_MULT = 0;
    NUM_OF_CURRENT_BOOT = 0;

}
void BootScheme::checkAndBoot(Ciphertext& cipher, bool condition, BootHelper& bootHelper, Parameter param) {
    if (condition) {
        PrintUtils::nprint("Run Boot in BootScheme::checkAndBoot -- before : " + to_string(cipher.logq), WANT_TO_PRINT);
        countBoot();
        // cout << "bootstrapping..." << endl;
        bootHelper.bootstrapping(cipher, param.logq, param.logQ, param.logT);
        PrintUtils::nprint("Run Boot in BootScheme::checkAndBoot -- after : " + to_string(cipher.logq), WANT_TO_PRINT);
    }
}
void BootScheme::multAndEqualWithBoot(Ciphertext& cipher1, Ciphertext& cipher2, BootHelper& bootHelper, Parameter param) {
    PrintUtils::nprint("multAndEqualWithBoot", WANT_TO_PRINT);
    if (cipher1.logq - param.logp < param.logq) {
        PrintUtils::nprint("before boot : " + to_string(cipher1.logq), WANT_TO_PRINT);
        PrintUtils::nprint("Bootstrapp Cipher1", WANT_TO_PRINT);
        bootHelper.bootstrapping(cipher1, param.logq, param.logQ, param.logT);
        PrintUtils::nprint("Bootstrapp Cipher2", WANT_TO_PRINT);
        bootHelper.bootstrapping(cipher2, param.logq, param.logQ, param.logT);
        PrintUtils::nprint("after boot : " + to_string(cipher1.logq), WANT_TO_PRINT);
    }
    countMult();
    multAndEqual(cipher1, cipher2);
}

void BootScheme::modDownToAndEqualModified(Ciphertext& cipher1, Ciphertext& cipher2, BootHelper& bootHelper, Parameter param) {
    checkAndBoot(cipher1, cipher1.logq < cipher2.logq, bootHelper, param);
    if (cipher1.logq < cipher2.logq) {
        modDownToAndEqual(cipher2, cipher1.logq);
    } else {
        modDownToAndEqual(cipher1, cipher2.logq);    
    }    
}

void BootScheme::squareAndEuqalWithBoot(Ciphertext& cipher, BootHelper& bootHelper, Parameter param) {
    countMult();
    checkAndBoot(cipher, cipher.logq - param.logp < param.logq, bootHelper, param);
    squareAndEqual(cipher);
}

void BootScheme::multByPolyAndEqualWithBoot(Ciphertext& cipher, ZZ* poly, BootHelper& bootHelper, Parameter param) {
    PrintUtils::nprint("in multByPolyAndEqualWithBoot : cipher.logq - param.logp = " + to_string(cipher.logq - param.logp) + " < param.logq = " + to_string(param.logq), WANT_TO_PRINT);
    checkAndBoot(cipher, cipher.logq - param.logp < param.logq, bootHelper, param);
    multByPolyAndEqual(cipher, poly, param.logp);
}

void BootScheme::showTotalCount() {
    cout << "Total Number of Multiplication = " << NUM_OF_MULT << endl;
    cout << "Total Number of Bootstrapping  = " << NUM_OF_BOOT << endl;
}

void BootScheme::showCurrentCount() {
    cout << "Number of Multiplication = " << NUM_OF_CURRENT_MULT << endl;
    cout << "Number of Bootstrapping  = " << NUM_OF_CURRENT_BOOT << endl;
}