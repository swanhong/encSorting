#include "test/TestPlain.h"
#include "test/TestEnc.h"
#include "test/TestBoot.h"
#include "test/TestSort.h"

int main() { 
    // ******************************
    // Parameters (long)
    // {logN, logQ, logp, logc, log2n, radix, logq, logT}
    // ******************************
    Parameter param = {11, 2000, 40, 40, 10, 32, 45, 4};
    Parameter bootTestParam = {7, 2000, 50, 50, 6, 8, 60, 4};
    Parameter sortingTestParamSmall = {9, 2000, 70, 70, 6, 64, 75, 4};
    Parameter sortingTestParam1 = {12, 1500, 30, 30, 10, 32, 35, 5};
    Parameter sortingTestParamBig = {15, 1200, 40, 40, 8, 16, 45, 5};
    Parameter sortingTestParamBig2 = {16, 1200, 40, 40, 15, 32, 45, 4};
    Parameter sortingTestParamBig3 = {15, 1200, 40, 40, 8, 16, 45, 5};
    Parameter sortingTestParam3 = {13, 2000, 30, 30, 12, 64, 35, 4};

    
    // ******************************
    // *** Test PlainSort
    // ******************************
    // TestPlain::showMasking(5, true);
    // TestPlain::showMasking(5, false);
    // TestPlain::showBitonicMergeMasking(2);
    // TestPlain::plainSort(5, false);
    // TestPlain::bitonicMerge(3, 5);

    // ******************************
    // *** Test algorithms for encrypted data
    // ******************************
    // TestEnc::approxSqrt(sortingTestParam1, 10);
    // TestEnc::minMax(sortingTestParam1, 10);
    // TestEnc::compAndSwap(sortingTestParamSmall, 15);

    // ******************************
    // *** Test algorithms for encrypted data with Bootstrapping
    // ******************************
    // TestBoot::approxSqrt(sortingTestParamSmall, 15);
    // TestBoot::minMax(sortingTestParam1, 12);
    // TestBoot::compAndSwap(sortingTestParamSmall, 10);
    // TestBoot::testSelfBitonicMerge(sortingTestParamSmall, 15);

    // ******************************
    // *** Check Parameters
    // ******************************
    // Parameter sortingTestParamSmall = {9, b , 70, 70, 6, 64, 75, 4};
    TestBoot::bootstrapping(param);
    // TestEnc::compAndSwap(sortingTestParamSmall, 10);
    // TestEnc::compAndSwap(sortingTestParamBig3, 15);
    // TestBoot::compAndSwap(sortingTestParamSmall, 15);
    // TestBoot::compAndSwapWithBoot(sortingTestParamSmall, 15);    

    // ******************************
    // *** Test EncSorting
    // ******************************
    // TestSort::sort(bootTestParam, 15);
    // TestSort::bitonicMerge(sortingTestParamSmall, 15);
    // TestSort::testMerge(sortingTestParamSmall, 15, 2);

    // {logN, logQ, logp, logc, log2n, radix, logq, logT}

    // Parameter param = {11, 2000, 50, 50, 10, 32, 60, 4};
    // srand(time(NULL));
	// SetNumThreads(8);
    
    // long n = 1 << param.log2n;
	
    // TimeUtils timeutils;
    // timeutils.start("KeyGen");
    // Ring ring(param.logN, param.logQ);
    // SecretKey secretKey(ring);
    // Scheme scheme(secretKey, ring);
    // scheme.addConjKey(secretKey);
    // scheme.addLeftRotKeys(secretKey);
    // scheme.addRightRotKeys(secretKey);
    // timeutils.stop("KeyGen");

    // BootHelper boothelper(param.log2n, param.radix, param.logc, scheme, ring, secretKey);

    // // //****** bootstrapping start

    // double* mvec = EvaluatorUtils::randomRealArray(n);
    
    // // for(long i = 0; i < n; i++) {
    // //     // mvec[i] =-0.25 + ((double)(i + 1) / (2*n + 1));
    // //     // cout << "mvec[" << i << "] = " << mvec[i] << endl;
    // //     // mvec[i] = 2 + (double) rand() / RAND_MAX / 100;
    // //     // mvec[i] += 3. + (double) (i % 3);
    // // }

	// Ciphertext cipher = scheme.encrypt(mvec, n, param.logp, param.logQ);

    // long logq = param.logq;
    // long logQ = param.logQ;
    // long logSlots = log2(cipher.n);
	// long logp = cipher.logp;

	// scheme.modDownToAndEqual(cipher, logq);
	// scheme.normalizeAndEqual(cipher);

	// TimeUtils time;

	// cipher.logq = logQ;
	// cipher.logp = logq;

	// for (long i = logSlots; i < ring.logNh; ++i) {
	// 	Ciphertext rot = scheme.leftRotateFast(cipher, (1 << i));
	// 	scheme.addAndEqual(cipher, rot);
	// }
	// scheme.divByPo2AndEqual(cipher, ring.logNh - logSlots);
	
	// Ciphertext part1, part2;

	// boothelper.coeffToSlot(part1, part2, cipher);

    // // complex<double>* dvec1 = scheme.decrypt(secretKey, part1);
    // // complex<double>* dvec2 = scheme.decrypt(secretKey, part2);
    
    // // boothelper.evalExpAndEqual(part1, part2, 4, 4, logq);
    // long logK = 6;
    // boothelper.evalSin2piAndEqual(part1, logK, logq);
	// boothelper.evalSin2piAndEqual(part2, logK, logq);
    // // complex<double>* dvec3 = scheme.decrypt(secretKey, part1);
    // // complex<double>* dvec4 = scheme.decrypt(secretKey, part2);

    // // for(long i = 0; i < n; i++) {
    // //     cout << dvec1[i].real() << " // " << dvec3[i].real() << " // " << endl;
    // //     cout << dvec2[i].real() << " // " << dvec4[i].real() << " // " << endl;
    // // }

    // boothelper.slotToCoeff(cipher, part1, part2);

    // cipher.logp = logp;

    // // //* ========================================================


	// // // evaluate x -> cos(2pi * x)
    // // double coeff[] = {1, -1, 2, 1};
    // // scheme.evalPolyAndEqual(cipher, param.logp, coeff, 0, 4);

    // // timeutils.start("eval cos");
    // // scheme.evalExpAndEqual(cipher, 4, 4);
	// // scheme.cos2piAndEqual(cipher, param.logp);
    // // timeutils.stop("eval cos");
    
    // // double* mvec = EvaluatorUtils::randomRealArray(n, 1./1024);
    // // for(long i = 0; i < n; i++) {
    // //     mvec[i] += 3. + (double) (i % 4);
    // // }
    // // Ciphertext cipher = scheme.encrypt(mvec, n, param.logp, param.logQ);
    
    // // cipher.logp = param.logq;
    // // Ciphertext cipher2 = cipher;
    // // // boothelper.evalExpAndEqual(cipher, cipher2, 4, 4, param.logq);
    // // boothelper.evalSin2piAndEqual(cipher, 8, param.logq);
    // // cipher.logp = param.logp;

	// // for(int i = 0; i < logK; i++) {
	// // 	scheme.squareAndEqual(cipher);
	// // 	scheme.multByConstAndEqual(cipher, 2.0, logp);
	// // 	scheme.reScaleByAndEqual(cipher, 2 * logp);
	// // 	scheme.addConstAndEqual(cipher, -1.0, logp);
	// // }	
    // // ************************************


    // // ************************************
    // // *** boothelper.evalExpAndEqual()
    // // ************************************
    // // long logT = 4;
    // // long logI = 4;

    // // cipher.logq = logQ;
	// // cipher.logp = logq;

    // // scheme.imultAndEqual(cipher); // i * theta
	// // cipher.logp += logI;
	// // scheme.divByPo2AndEqual(cipher, logT); // i * theta / (2^(logT + logI))
	// // scheme.exp2piAndEqual(cipher, logq + logI); // exp(2pi*i * theta / (2^(logT + logI)))
	// // for (long i = 0; i < logI + logT; ++i) {
	// // 	scheme.squareAndEqual(cipher);
	// // 	scheme.reScaleByAndEqual(cipher, logq + logI);
	// // } // exp(2pi * i * theta) = cos(2pi * theta) + sin(2pi * theta) * i
    // // // ************************************
    // // Ciphertext tmp1 = scheme.conjugate(cipher);
	// // scheme.subAndEqual(cipher, tmp1);
	// // scheme.imultAndEqual(cipher);
	// // scheme.negateAndEqual(cipher); // 2 * sin(2pi * theta)
	// // RR c = 0.25 / to_RR(M_PI);
	// // scheme.multByConstAndEqual(cipher, c, logq + logI); // 1/2pi * (sin(2pi * theta))
	// // scheme.reScaleByAndEqual(cipher, logq + 2 * logI);

    // complex<double>* dvec = scheme.decrypt(secretKey, cipher);
    // // complex<double>* dvec1 = scheme.decrypt(secretKey, part1);
    // // complex<double>* dvec2 = scheme.decrypt(secretKey, part2);
    // // PrintUtils::printArrays(mvec, dvec, n);
    // cout << "consumed logQ = " << logQ - cipher.logq << endl;
    // PrintUtils::averageDifference(mvec, dvec, n);

    return 0;
}