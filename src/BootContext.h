#ifndef BOOTCONTEXT_H_
#define BOOTCONTEXT_H_

#include "../HEAAN/src/HEAAN.h"
#include "../bootsrc/new_bootstrapping.h"
#include "Parameter.h"
#include "SortingAlgorithm.h"
#include "EncContext.h"

void multWithBootAndEqual(Scheme& scheme, Ciphertext& cipher1, Ciphertext& cipher2, BootHelper boothelper, Parameter parameter);

void fcnEncSqrtWithBoot(Ciphertext& outCipher, const Ciphertext& inCipher, Parameter parameter, long d, Scheme& scheme, BootHelper& boothelper);

void fcnEncMaxMinWithBoot(Ciphertext& maxCipher, Ciphertext& minCipher, Ciphertext& input1, Ciphertext& input2, Parameter parameter, long iter, Scheme& scheme, BootHelper& boothelper);

void fcnEncCompAndSwapWithBoot(Ciphertext& outCipher, const Ciphertext& inCipher, double* mask, long dist, Parameter parameter, long iter, Scheme& scheme, BootHelper& boothelper);

long encBatcherOddEvenSortRec(Ciphertext& sortedCipher, const Ciphertext& inCipher, long logNum, long logJump, long loc, Parameter parameter, long iter, double** mask, Scheme& scheme, BootHelper& boothelper);

void encBatcherOddEvenSort(Ciphertext& sortedCipher, const Ciphertext& inCipher, Parameter parameter, long iter, Scheme& scheme, BootHelper& boothelper);

void encBatcherOddEvenSortWithDecrypt(Ciphertext& sortedCipher, const Ciphertext& inCipher, Parameter parameter, long iter, Scheme& scheme, BootHelper& boothelper, SecretKey& secretKey);

long encBatcherOddEvenSortRecWithDecrypt(Ciphertext& sortedCipher, const Ciphertext& inCipher, long logNum, long logJump, long loc, Parameter parameter, long iter, double** mask, Scheme& scheme, BootHelper& boothelper, SecretKey& secretKey);

void fcnEncCompAndSwapWithBootAndDecrypt(Ciphertext& outCipher, Ciphertext& inCipher, double* mask, long dist, Parameter parameter, long iter, Scheme& scheme, BootHelper& boothelper, SecretKey& secretKey);

void fcnEncMaxMinWithBootAndDecrypt(Ciphertext& maxCipher, Ciphertext& minCipher, Ciphertext& input1, Ciphertext& input2, Parameter parameter, long iter, Scheme& scheme, BootHelper& boothelper, SecretKey& secretKey);

void fcnEncSqrtWithBootWithDecrypt(Ciphertext& outCipher, Ciphertext& inCipher, Parameter parameter, long d, Scheme& scheme, BootHelper& boothelper, SecretKey& secretKey);
// ======================================================================
// ======================================================================
// ======================================================================
// ======================================================================
// ======================================================================


#endif // !BOOTCONTEXT_H_