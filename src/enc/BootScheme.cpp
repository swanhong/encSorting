#include "BootScheme.h"


void BootScheme::checkAndBoot(Ciphertext& cipher, bool condition, BootHelper& bootHelper, Parameter param) {
    if (condition) {
        cout << "Ron Boot -- before : " << cipher.logq << endl;
        PrintUtils::nprint("Bootstrapp Cipher", WANT_TO_PRINT);
        bootHelper.bootstrapping(cipher, param.logq, param.logQ, param.logT);
        cout << "run boot -- after : " << cipher.logq << endl;
    }
}
void BootScheme::multAndEqualWithBoot(Ciphertext& cipher1, Ciphertext& cipher2, BootHelper& bootHelper, Parameter param) {
    cout << "multAndEqualWithBoot" << endl;
    if (cipher1.logq - param.logp < param.logq) {
        cout << "before boot : " << cipher1.logq << endl;
        PrintUtils::nprint("Bootstrapp Cipher1", WANT_TO_PRINT);
        bootHelper.bootstrapping(cipher1, param.logq, param.logQ, param.logT);
        PrintUtils::nprint("Bootstrapp Cipher2", WANT_TO_PRINT);
        bootHelper.bootstrapping(cipher2, param.logq, param.logQ, param.logT);
        cout << "after boot : " << cipher1.logq << endl;
    }
    PrintUtils::nprint("Mult", WANT_TO_PRINT);
    multAndEqual(cipher1, cipher2);
}

void BootScheme::modDownToAndEqualModified(Ciphertext& cipher1, Ciphertext& cipher2, BootHelper& bootHelper, Parameter param) {
    cout << "before checkboot : " << "cipher1.logq = " << cipher1.logq << " vs cipher2.logq = " << cipher2.logq << endl;
    checkAndBoot(cipher1, cipher1.logq < cipher2.logq, bootHelper, param);
    if (cipher1.logq < cipher2.logq) {
        modDownToAndEqual(cipher2, cipher1.logq);
    } else {
        modDownToAndEqual(cipher1, cipher2.logq);    
    }    
}

void BootScheme::squareAndEuqalWithBoot(Ciphertext& cipher, BootHelper& bootHelper, Parameter param) {
    checkAndBoot(cipher, cipher.logq - param.logp < param.logq, bootHelper, param);
    squareAndEqual(cipher);
}

void BootScheme::multByPolyAndEqualWithBoot(Ciphertext& cipher, ZZ* poly, BootHelper& bootHelper, Parameter param) {
    cout << "in multByPolyAndEqualWithBoot : cipher.logq - param.logp = " << cipher.logq - param.logp << " < param.logq = " << param.logq << endl;
    checkAndBoot(cipher, cipher.logq - param.logp < param.logq, bootHelper, param);
    multByPolyAndEqual(cipher, poly, param.logp);
}

