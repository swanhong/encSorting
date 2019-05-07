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

void BootScheme::checkLevelAndBoot(Ciphertext& cipher, long level, BootHelper& bootHelper, Parameter param) {
    bool condition = cipher.logq - level * param.logp < param.logq;
    checkAndBoot(cipher, condition, bootHelper, param);
}

void BootScheme::checkModulusAndBoot(Ciphertext& cipher, long mod, BootHelper& bootHelper, Parameter param) {
    bool condition = cipher.logq - mod < param.logq;
    checkAndBoot(cipher, condition, bootHelper, param);
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
    if (cipher1.logq != cipher2.logq) {
        cout << " =========== error, mult for different logq in BootScheme::multAndEqualWithBoot ============" << endl;
    }
    
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

Ciphertext BootScheme::multByVectorWithBoot(Ciphertext& cipher, double* mask, long loga, Ring& ring, BootHelper& bootHelper, Parameter param) {
    Ciphertext res = cipher;
    multByVectorAndEqualWithBoot(res, mask, loga, ring, bootHelper, param);
    return res;
}

void BootScheme::multByVectorAndEqualWithBoot(Ciphertext& cipher, double* mask, long loga, Ring& ring, BootHelper& bootHelper, Parameter param) {
    ZZ* maskPoly = new ZZ[1 << param.logN];
    ring.encode(maskPoly, mask, cipher.n, param.logp);
    multByPolyAndEqualWithBoot(cipher, maskPoly, loga, bootHelper, param);
}

Ciphertext BootScheme::multByPolyWithBoot(Ciphertext& cipher, ZZ* poly, long loga, BootHelper& bootHelper, Parameter param) {
    Ciphertext res = cipher;
    multByPolyAndEqualWithBoot(res, poly, loga, bootHelper, param);
    return res;
}

void BootScheme::multByPolyAndEqualWithBoot(Ciphertext& cipher, ZZ* poly, long loga, BootHelper& bootHelper, Parameter param) {
    PrintUtils::nprint("in multByPolyAndEqualWithBoot : cipher.logq - param.logp = " + to_string(cipher.logq - param.logp) + " < param.logq = " + to_string(param.logq), WANT_TO_PRINT);
    checkModulusAndBoot(cipher, loga, bootHelper, param);
    // multByPolyAndEqual(cipher, poly, param.logp);
    multByPolyAndEqual(cipher, poly, loga);
}

Ciphertext BootScheme::leftRotateConditional(Ciphertext& cipher, long r, bool condition) {
    Ciphertext res = cipher;
    leftRotateAndEqualConditional(res, r, condition);
    return res;
}

Ciphertext BootScheme::rightRotateConditional(Ciphertext& cipher, long r, bool condition) {
    Ciphertext res = cipher;
    rightRotateAndEqualConditional(res, r, condition);
    return res;
}

void BootScheme::leftRotateAndEqualConditional(Ciphertext& cipher, long r, bool condition) {
    if (condition) leftRotateFastAndEqual(cipher, r);
    else rightRotateFastAndEqual(cipher, r);
}

void BootScheme::rightRotateAndEqualConditional(Ciphertext& cipher, long r, bool condition) {
    if (condition) rightRotateFastAndEqual(cipher, r);
    else leftRotateFastAndEqual(cipher, r);
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