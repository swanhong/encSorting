/*
* Copyright (c) by CryptoLab inc.
* This program is licensed under a
* Creative Commons Attribution-NonCommercial 3.0 Unported License.
* You should have received a copy of the license along with this
* work.  If not, see <http://creativecommons.org/licenses/by-nc/3.0/>.
*/
#ifndef HEAAN_SERIALIZATIONUTILS_H_
#define HEAAN_SERIALIZATIONUTILS_H_

#include <iostream>

#include "Ciphertext.h"
#include "Plaintext.h"
#include "Ring.h"
#include "Scheme.h"
#include "SecretKey.h"
#include "Key.h"

using namespace std;
using namespace NTL;

class SerializationUtils {
public:

	static void writeCiphertext(Ciphertext& ciphertext, string path);
	static Ciphertext readCiphertext(string path);

	static void writeKey(Key* key, string path);
	static Key* readKey(string path);

	static void writeSecretKey(SecretKey& sk, string path);
	static SecretKey readSecretKey(string path);
};

#endif /* SERIALIZATIONUTILS_H_ */
