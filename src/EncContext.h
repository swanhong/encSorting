#ifndef ENCCONTEXT_H_
#define ENCCONTEXT_H_

#include "../HEAAN/src/HEAAN.h"

void fcnDecryptAndPrint(string str, Ciphertext& cipher, Scheme& scheme, SecretKey& secretKey);
void fcnDecryptAndPrintTwo(string str, Ciphertext& cipher1, Ciphertext& cipher2, Scheme& scheme, SecretKey& secretKey);

void fcnEncComputeSqrt(Ciphertext& outCipher, const Ciphertext& inCipher, long logp, long d, Scheme scheme);

void fcnEncMaxMin(Ciphertext& maxCipher, Ciphertext& minCipher, Ciphertext& input1, Ciphertext& input2, long logp, long iter, Scheme scheme);

void fcnGenEncAllMasking(long log2n, double**& mask);

long fcnGenEncMaskingRec(long log2n, long logNum, long logJump, double**& mask, long loc);

double* fcnGenEncMaskingComp(long length, long jump);

double* fcnGenEncMaskingMerge(long length, long num, long jump);

void fcnEncCompAndSwap(Ciphertext& outCipher, const Ciphertext& inCipher, double* mask, long dist, long logp, long iter, Scheme scheme);

#endif // !ENCCONTEXT_H_