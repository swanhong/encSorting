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
        // bootHelper.bootstrapping(cipher, param.logq, param.logQ, param.logT);
        cout << "bootstrapping...." << endl;
        bootHelper.bootstrapping_cos(cipher, param.logq, param.logQ, 5);
        PrintUtils::nprint("Run Boot in BootScheme::checkAndBoot -- after : " + to_string(cipher.logq), WANT_TO_PRINT);
    }
}

Ciphertext BootScheme::multWithBoot(Ciphertext& cipher1, Ciphertext& cipher2, BootHelper& bootHelper, Parameter param) {
    Ciphertext res = cipher1;
    multAndEqualWithBoot(res, cipher2, bootHelper, param);
    return res;
}

void BootScheme::multAndEqualWithBoot(Ciphertext& cipher1, Ciphertext& cipher2, BootHelper& bootHelper, Parameter param) {
    checkAndBoot(cipher1, cipher1.logq - param.logp < param.logq,bootHelper,param);
    checkAndBoot(cipher2, cipher2.logq - param.logp < param.logq,bootHelper,param);
    countMult();
    multAndEqual(cipher1, cipher2);
}

void BootScheme::modDownToAndEqualModified(Ciphertext& cipher1, Ciphertext& cipher2, BootHelper& bootHelper, Parameter param) {
    if (cipher1.logq < cipher2.logq) {
        modDownToAndEqual(cipher2, cipher1.logq);
    } else {
        modDownToAndEqual(cipher1, cipher2.logq);    
    }    
}

Ciphertext BootScheme::squareWithBoot(Ciphertext& cipher, BootHelper& bootHelper, Parameter param) {
    Ciphertext res = cipher;
    squareAndEuqalWithBoot(res, bootHelper, param);
    return res;
}
void BootScheme::squareAndEuqalWithBoot(Ciphertext& cipher, BootHelper& bootHelper, Parameter param) {
    countMult();
    checkAndBoot(cipher, cipher.logq - param.logp < param.logq, bootHelper, param);
    squareAndEqual(cipher);
}

Ciphertext BootScheme::multByPolyWithBoot(Ciphertext& cipher, ZZ* poly, BootHelper& bootHelper, Parameter param) {
    Ciphertext res = cipher;
    multByPolyAndEqualWithBoot(res, poly, bootHelper, param);
    return res;
}

void BootScheme::multByPolyAndEqualWithBoot(Ciphertext& cipher, ZZ* poly, BootHelper& bootHelper, Parameter param) {
    PrintUtils::nprint("in multByPolyAndEqualWithBoot : cipher.logq - param.logp = " + to_string(cipher.logq - param.logp) + " < param.logq = " + to_string(param.logq), WANT_TO_PRINT);
    checkAndBoot(cipher, cipher.logq - param.logp < param.logq, bootHelper, param);
    multByPolyAndEqual(cipher, poly, param.logp);
}

void BootScheme::nomalizeAndEuqal(Ciphertext& cipher) {
    long dup = log2(cipher.N / cipher.n / 2);
    // cout << "nomalize with dup = " << dup << endl;
    for (int i = 0; i < dup; i++) {
        // cout << "rot by " << cipher.n * (1 << i) << endl;
        Ciphertext rot = leftRotateFast(cipher, cipher.n * (1 << i));
        addAndEqual(cipher, rot);
    }
    // cout << "div by " << dup << endl;
    divByPo2AndEqual(cipher, dup);
    
}

void BootScheme::decryptAndPrint(std::string str, SecretKey& secretKey, Ciphertext& cipher) {
    complex<double>* dvec = decrypt(secretKey, cipher);
    PrintUtils::printSingleArray(str, dvec, cipher.n);
}

void BootScheme::showTotalCount() {
    cout << "Total Number of Multiplication = " << NUM_OF_MULT << endl;
    cout << "Total Number of Bootstrapping  = " << NUM_OF_BOOT << endl;
}

void BootScheme::showCurrentCount() {
    cout << "Number of Multiplication = " << NUM_OF_CURRENT_MULT << endl;
    cout << "Number of Bootstrapping  = " << NUM_OF_CURRENT_BOOT << endl;
}