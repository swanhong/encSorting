// #include "../bootsrc/new_bootstrapping.h"
#include "Parameter.h"
#include "BootContext.h"
#include "stdlib.h"


class TestBootstrapping {
public:
    static void bootstrapping_test(Parameter parameter);

	static void bootstrapping_test_with_mult(Parameter parameter, int iter);
	
	static void testSqrtWithBoot(Parameter parameter, long iter);
	
	static void testMaxMinWithBoot(Parameter parameter, long iter);

	static void testEncCompAndSwapWithBoot(Parameter parameter, long iter);

	static void testEncSort(Parameter parameter, long iter);

	static void testEncSortWithDecrypt(Parameter parameter, long iter);

	static void testMaxMinWithBootAndDecrypt(Parameter parameter, long iter);

	static void testSqrtWithBootAndDecrypt(Parameter parameter, long iter);

	static void testEncCompAndSwapWithBootAndDecrypt(Parameter parameter, long iter);
};

void TestBootstrapping::bootstrapping_test(Parameter parameter) {
	// HE parameter //
	long logN = parameter.logN;
	long logQ = parameter.logQ;
	long logp = parameter.logp;
	long logc = parameter.logc;

	// Decomposition related parameter //
	long log2n = parameter.log2n;
	long radix = parameter.radix;

	// Bootstrapping parameter //
	long logq = parameter.logq;
	long logT = parameter.logT;

	long n = 1 << log2n;

	cout << "\n***************************" << endl;
	cout << "Test for Improved Bootstrapping" << endl;
	cout << "logN = " << logN << ", logQ = " << logQ << ", logp = " << logp << ", logc = " << logc << endl;
	cout << "slots = " << n << ", radix = " << radix << ", logq = " << logq << ", logT = " << logT << endl;
	cout << "***************************" << endl;
	cout << endl;

	TimeUtils timeutils;
	timeutils.start("KeyGen");
	Ring ring(logN, logQ);
	SecretKey secretKey(ring);
	Scheme scheme(secretKey, ring);
	scheme.addConjKey(secretKey);
	scheme.addLeftRotKeys(secretKey);
	timeutils.stop("KeyGen");

	timeutils.start("Bootstrapping Helper construct");
	BootHelper boothelper(log2n, radix, logc, scheme, ring, secretKey);
	timeutils.stop("Bootstrapping Helper construct");

	complex<double>* mvec = EvaluatorUtils::randomComplexArray(n);
	// complex<double>* mvec = new complex<double>[n];
	// for(int i = 0; i < n; i++) {
	// 	mvec[i] *= 0.0000001;
	// }
	mvec[0] = complex<double>(0.5, 0.3);
	
	Ciphertext cipher = scheme.encrypt(mvec, n, logp, logQ);
	mvec = scheme.decrypt(secretKey, cipher);

	complex<double>* dvec;
	
	// scheme.multByConstAndEqual(cipher, 0.01, logp);
	// scheme.reScaleByAndEqual(cipher, logp);
	
	// Bootstrapping //
	timeutils.start("Improved bootstrapping");
	boothelper.bootstrapping(cipher, logq, logQ, logT);
	timeutils.stop("Improved bootstrapping");

	// scheme.multByConstAndEqual(cipher, 100, logp);
	// scheme.reScaleByAndEqual(cipher, logp);

	cout << "* Befor logQ = " << logq << endl;
	cout << "* After logQ = " << cipher.logq << endl;

	// Print Result and Difference //
	dvec = scheme.decrypt(secretKey, cipher);
	cout << "log2(avg of error) = " << diff(mvec, dvec, n) << endl;

	// for(int i = 0; i < n; i++)
	// {
	// 	cout << i << " : " << mvec[i] << " ~~ " << dvec[i] << endl;
	// }
	
	if(mvec != NULL) delete[] mvec;
	if(dvec != NULL) delete[] dvec;
	return;
}

void TestBootstrapping::bootstrapping_test_with_mult(Parameter parameter, int iter) {
	// HE parameter //
	long logN = parameter.logN;
	long logQ = parameter.logQ;
	long logp = parameter.logp;
	long logc = parameter.logc;

	// Decomposition related parameter //
	long log2n = parameter.log2n;
	long radix = parameter.radix;

	// Bootstrapping parameter //
	long logq = parameter.logq;
	long logT = parameter.logT;

	long n = 1 << log2n;

	cout << "\n***************************" << endl;
	cout << "Test for Improved Bootstrapping with mult" << endl;
	cout << "logN = " << logN << ", logQ = " << logQ << ", logp = " << logp << ", logc = " << logc << endl;
	cout << "slots = " << n << ", radix = " << radix << ", logq = " << logq << ", logT = " << logT << endl;
	cout << "***************************" << endl;
	cout << endl;

	TimeUtils timeutils;
	timeutils.start("KeyGen");
	Ring ring(logN, logQ);
	SecretKey secretKey(ring);
	Scheme scheme(secretKey, ring);
	scheme.addConjKey(secretKey);
	scheme.addLeftRotKeys(secretKey);
	timeutils.stop("KeyGen");

	timeutils.start("Bootstrapping Helper construct");
	BootHelper boothelper(log2n, radix, logc, scheme, ring, secretKey);
	timeutils.stop("Bootstrapping Helper construct");

	complex<double>* mvec1 = EvaluatorUtils::randomComplexArray(n, 10);
	
	Ciphertext cipher1 = scheme.encrypt(mvec1, n, logp, logQ);

	// complex<double>* mvec2 = EvaluatorUtils::randomComplexArray(n, 0.7);
	// Ciphertext cipher2 = scheme.encrypt(mvec2, n, logp, logQ);
	
	complex<double>* mvec2 = new complex<double>[n];
	for (int i = 0; i < n; i++) {
		mvec2[i] = 1;
	}
	Ciphertext cipher2 = scheme.encrypt(mvec2, n, logp, logQ);

	cout << "======= initial values ======= " << endl;
	for(int i = 0; i < n; i++) {
		cout << i << " : " << mvec1[i] << " ~~ " << mvec2[i] << endl;
	}
	
	
	complex<double>* dvec;
	
	// mult mult mult //

	for(int i = 0; i < iter; i++) {
		cout << cipher1.logq - logp << " vs " << logq << endl;
		if (cipher1.logq - logp <  logq) {
		// if (1) {
			cout << "Bootstrapping" << endl;
			cout << "...before bootstrapping : " << endl;
			cout << "logq = " << cipher1.logq << endl;
			complex<double>* dvec = scheme.decrypt(secretKey, cipher1);

			// scheme.multByConstAndEqual(cipher1, 1000000.0, logp);
			// scheme.multByConstAndEqual(cipher2, 1000000.0, logp);
			// scheme.reScaleByAndEqual(cipher1, logp);
			// scheme.reScaleByAndEqual(cipher2, logp);

			boothelper.bootstrapping(cipher1, logq, logQ, logT);
			boothelper.bootstrapping(cipher2, logq, logQ, logT);

			// scheme.multByConstAndEqual(cipher1, 0.000001, logp);
			// scheme.multByConstAndEqual(cipher2, 0.000001, logp);
			// scheme.reScaleByAndEqual(cipher1, logp);
			// scheme.reScaleByAndEqual(cipher2, logp);

			complex<double>* dvec2 = scheme.decrypt(secretKey, cipher1);
			cout << "...after bootstrapping : " << endl;
			cout << "logq = " << cipher1.logq << endl;
			for(int i = 0; i < cipher1.n; i++) {
				if (cipher1.n > 100) {
					if (i % 100 == 0) {
						cout << i << " : " << dvec[i] << " ~~ " << dvec2[i] << endl;
					}        
				} else {
					cout << i << " : " << dvec[i] << " ~~ " << dvec2[i] << endl;
				}        
			}
		}
		cout << i << "th mult" << endl;
		scheme.multAndEqual(cipher1, cipher2);
		scheme.reScaleByAndEqual(cipher1, logp);
		scheme.modDownByAndEqual(cipher2, logp);
		for(int j = 0; j < n; j++) {
			mvec1[j] *= mvec2[j];
		}
		// scheme.modDownByAndEqual(cipher1, logp);
		// scheme.modDownByAndEqual(cipher2, logp);

		complex<double>* dvec = scheme.decrypt(secretKey, cipher1);
		for(int i = 0; i < cipher1.n; i++) {
			if (cipher1.n > 100) {
				if (i % 100 == 0) {
					cout << i << " : " << mvec1[i] << " ~~ " << dvec[i] << endl;
				}        
			} else {
				cout << i << " : " << dvec[i] << endl;
			}        
		}
	}

	// fcnDecryptAndPrint("cipher : ", cipher1, scheme, secretKey);

	// Bootstrapping //
	// timeutils.start("Improved bootstrapping");
	// boothelper.bootstrapping(cipher, logq, logQ, logT);
	// timeutils.stop("Improved bootstrapping");

	// cout << "* Befor logQ = " << logq << endl;
	// cout << "* After logQ = " << cipher.logq << endl;

	// Print Result and Difference //
	dvec = scheme.decrypt(secretKey, cipher1);
	// cout << "log2(avg of error) = " << diff(mvec, dvec, n) << endl;

	cout << " ========== final output =========== " << endl;
	for(int i = 0; i < n; i++)
	{
		cout << i << " : " << mvec1[i] << " ~simeq~ " << dvec[i] << endl;
	}
	
	// if(mvec != NULL) delete[] mvec;
	// if(dvec != NULL) delete[] dvec;
	return;
}


void TestBootstrapping::testSqrtWithBoot(Parameter parameter, long iter) {
	// HE parameter //
	long logN = parameter.logN;
	long logQ = parameter.logQ;
	long logp = parameter.logp;
	long logc = parameter.logc;

	// Decomposition related parameter //
	long log2n = parameter.log2n;
	long radix = parameter.radix;

	// Bootstrapping parameter //
	long logq = parameter.logq;
	long logT = parameter.logT;

	long n = 1 << log2n;

	cout << "\n***************************" << endl;
	cout << "Test for Sqrt with Boot" << endl;
	cout << "logN = " << logN << ", logQ = " << logQ << ", logp = " << logp << ", logc = " << logc << endl;
	cout << "slots = " << n << ", radix = " << radix << ", logq = " << logq << ", logT = " << logT << endl;
	cout << "***************************" << endl;
	cout << endl;

	TimeUtils timeutils;
	timeutils.start("KeyGen");
	Ring ring(logN, logQ);
	SecretKey secretKey(ring);
	Scheme scheme(secretKey, ring);
	scheme.addConjKey(secretKey);
	scheme.addLeftRotKeys(secretKey);
	timeutils.stop("KeyGen");

	timeutils.start("Bootstrapping Helper construct");
	BootHelper boothelper(log2n, radix, logc, scheme, ring, secretKey);
	timeutils.stop("Bootstrapping Helper construct");

	double* mvec = new double[n];
    for(int i = 0; i < n; i++) {
        mvec[i] = 0.00001 + (double) i / n;
    }
    
	Ciphertext cipher = scheme.encrypt(mvec, n, logp, logQ);

    Ciphertext sqrtCipher;
    timeutils.start("Sqrt");
    // fcnEncComputeSqrt(sqrtCipher, cipher, logp, 5, scheme);
	fcnEncSqrtWithBoot(sqrtCipher, cipher, parameter, iter, scheme, boothelper);
    timeutils.stop("Sqrt");

    // Print Result and Difference //	
	complex<double>* dvec = scheme.decrypt(secretKey, sqrtCipher);
    for(int i = 0; i < n; i++) {
		double ans = dvec[i].real();
        cout << i << " : " << mvec[i] << " ~simeq~ "
			<< ans * ans << " = " << ans << " ^2" << endl;
    }
}

void TestBootstrapping::testMaxMinWithBoot(Parameter parameter, long iter) {
	// HE parameter //
	long logN = parameter.logN;
	long logQ = parameter.logQ;
	long logp = parameter.logp;
	long logc = parameter.logc;

	// Decomposition related parameter //
	long log2n = parameter.log2n;
	long radix = parameter.radix;

	// Bootstrapping parameter //
	long logq = parameter.logq;
	long logT = parameter.logT;

	long n = 1 << log2n;

	cout << "\n***************************" << endl;
	cout << "Test for MaxMin with Boot" << endl;
	cout << "logN = " << logN << ", logQ = " << logQ << ", logp = " << logp << ", logc = " << logc << endl;
	cout << "slots = " << n << ", radix = " << radix << ", logq = " << logq << ", logT = " << logT << endl;
	cout << "***************************" << endl;
	cout << endl;

	TimeUtils timeutils;
	timeutils.start("KeyGen");
	Ring ring(logN, logQ);
	SecretKey secretKey(ring);
	Scheme scheme(secretKey, ring);
	scheme.addConjKey(secretKey);
	scheme.addLeftRotKeys(secretKey);
	timeutils.stop("KeyGen");

	timeutils.start("Bootstrapping Helper construct");
	BootHelper boothelper(log2n, radix, logc, scheme, ring, secretKey);
	timeutils.stop("Bootstrapping Helper construct");

	double* mvec1 = new double[n];
    for(int i = 0; i < n; i++) {
        mvec1[i] = 0.00001 + (double) i / n;
    }

    double* mvec2 = new double[n];
    for(int i = 0; i < n; i++) {
        mvec2[i] = 1.0 - (double) i / n;
    }
    
	Ciphertext cipher1 = scheme.encrypt(mvec1, n, parameter.logp, parameter.logQ);
    Ciphertext cipher2 = scheme.encrypt(mvec2, n, parameter.logp, parameter.logQ);
    
    Ciphertext minCipher, maxCipher;
	timeutils.start("MaxMin with Bootstrapping");
    fcnEncMaxMinWithBoot(maxCipher, minCipher, cipher1, cipher2, parameter, iter, scheme, boothelper);
    timeutils.stop("MaxMin with Bootstrapping");

	timeutils.start("Decryption");
    complex<double>* dvec1 = scheme.decrypt(secretKey, maxCipher);
    complex<double>* dvec2 = scheme.decrypt(secretKey, minCipher);
    // StringUtils::compare(mvec, dmult, n, "mult");
	timeutils.stop("Decryption");

    for(int i = 0; i < n; i++) {
        cout << i << " : " << 
            mvec1[i] << ", " << mvec2[i] << " /// -> /// " <<
            dvec1[i].real() << ", " << dvec2[i].real() << endl;
    }
	
	return;
}

void TestBootstrapping::testMaxMinWithBootAndDecrypt(Parameter parameter, long iter) {
	// HE parameter //
	long logN = parameter.logN;
	long logQ = parameter.logQ;
	long logp = parameter.logp;
	long logc = parameter.logc;

	// Decomposition related parameter //
	long log2n = parameter.log2n;
	long radix = parameter.radix;

	// Bootstrapping parameter //
	long logq = parameter.logq;
	long logT = parameter.logT;

	long n = 1 << log2n;

	cout << "\n***************************" << endl;
	cout << "Test for MaxMin with Boot" << endl;
	cout << "logN = " << logN << ", logQ = " << logQ << ", logp = " << logp << ", logc = " << logc << endl;
	cout << "slots = " << n << ", radix = " << radix << ", logq = " << logq << ", logT = " << logT << endl;
	cout << "***************************" << endl;
	cout << endl;

	TimeUtils timeutils;
	timeutils.start("KeyGen");
	Ring ring(logN, logQ);
	SecretKey secretKey(ring);
	Scheme scheme(secretKey, ring);
	scheme.addConjKey(secretKey);
	scheme.addLeftRotKeys(secretKey);
	timeutils.stop("KeyGen");

	timeutils.start("Bootstrapping Helper construct");
	BootHelper boothelper(log2n, radix, logc, scheme, ring, secretKey);
	timeutils.stop("Bootstrapping Helper construct");

	double* mvec1 = new double[n];
	double* mvec2 = new double[n];
    for(int i = 0; i < n; i++) {
        mvec1[i] = (double) rand() / RAND_MAX;
		mvec2[i] = (double) rand() / RAND_MAX;
    }
    
	Ciphertext cipher1 = scheme.encrypt(mvec1, n, parameter.logp, parameter.logQ);
    Ciphertext cipher2 = scheme.encrypt(mvec2, n, parameter.logp, parameter.logQ);
    
    Ciphertext minCipher, maxCipher;
	timeutils.start("MaxMin with Bootstrapping");
    fcnEncMaxMinWithBootAndDecrypt(maxCipher, minCipher, cipher1, cipher2, parameter, iter, scheme, boothelper, secretKey);
    timeutils.stop("MaxMin with Bootstrapping");

	timeutils.start("Decryption");
    complex<double>* dvec1 = scheme.decrypt(secretKey, maxCipher);
    complex<double>* dvec2 = scheme.decrypt(secretKey, minCipher);
    // StringUtils::compare(mvec, dmult, n, "mult");
	timeutils.stop("Decryption");

    for(int i = 0; i < n; i++) {
        cout << i << " : " << 
            mvec1[i] << ", " << mvec2[i] << " /// -> /// " <<
            dvec1[i].real() << ", " << dvec2[i].real() << endl;
    }
	
	return;
}

void TestBootstrapping::testEncCompAndSwapWithBoot(Parameter parameter, long iter) {
         // HE parameter //
	long logN = parameter.logN;
	long logQ = parameter.logQ;
	long logp = parameter.logp;
	long logc = parameter.logc;

	// Decomposition related parameter //
	long log2n = parameter.log2n;
	long radix = parameter.radix;

	// Bootstrapping parameter //
	long logq = parameter.logq;
	long logT = parameter.logT;

	long n = 1 << log2n;

	cout << "\n***************************" << endl;
	cout << "Test for Comparison And Swap with Boot" << endl;
	cout << "logN = " << logN << ", logQ = " << logQ << ", logp = " << logp << ", logc = " << logc << endl;
	cout << "slots = " << n << ", radix = " << radix << ", logq = " << logq << ", logT = " << logT << endl;
	cout << "***************************" << endl;
	cout << endl;

	TimeUtils timeutils;
	timeutils.start("KeyGen");
	Ring ring(logN, logQ);
	SecretKey secretKey(ring);
	Scheme scheme(secretKey, ring);
	scheme.addConjKey(secretKey);
	scheme.addLeftRotKeys(secretKey);
	scheme.addRightRotKeys(secretKey);
	timeutils.stop("KeyGen");

	timeutils.start("Bootstrapping Helper construct");
	BootHelper boothelper(log2n, radix, logc, scheme, ring, secretKey);
	timeutils.stop("Bootstrapping Helper construct");

    double* mvec = new double[n];
    for(int i = 0; i < n; i++) {
        // mvec[n - 1 - i] = 0.00001 + (double) i / n;
		mvec[i] = (double) rand() / RAND_MAX;
    }
    
	Ciphertext cipher = scheme.encrypt(mvec, n, logp, logQ);

	Ciphertext cipher2;

    double* mask = genMaskingComp(n, 1);
    
    timeutils.start("CompAndSwap");
    fcnEncCompAndSwapWithBoot(cipher2, cipher, mask, 1, parameter, iter, scheme, boothelper);
    timeutils.stop("CompAndSwap");

	complex<double>* dvec = scheme.decrypt(secretKey, cipher2);
    // StringUtils::compare(mvec, dmult, n, "mult");
    for(int i = 0; i < n; i++) {
        cout << i << " : " << mvec[i] << ", " << dvec[i].real() << endl;
    }
}

void TestBootstrapping::testEncSort(Parameter parameter, long iter) {
         // HE parameter //
	long logN = parameter.logN;
	long logQ = parameter.logQ;
	long logp = parameter.logp;
	long logc = parameter.logc;

	// Decomposition related parameter //
	long log2n = parameter.log2n;
	long radix = parameter.radix;

	// Bootstrapping parameter //
	long logq = parameter.logq;
	long logT = parameter.logT;

	long n = 1 << log2n;

	cout << "\n***************************" << endl;
	cout << "Test for EncSort" << endl;
	cout << "logN = " << logN << ", logQ = " << logQ << ", logp = " << logp << ", logc = " << logc << endl;
	cout << "slots = " << n << ", radix = " << radix << ", logq = " << logq << ", logT = " << logT << endl;
	cout << "***************************" << endl;
	cout << endl;

	TimeUtils timeutils;
	timeutils.start("KeyGen");
	Ring ring(logN, logQ);
	SecretKey secretKey(ring);
	Scheme scheme(secretKey, ring);
	scheme.addConjKey(secretKey);
	scheme.addLeftRotKeys(secretKey);
	scheme.addRightRotKeys(secretKey);
	timeutils.stop("KeyGen");

	timeutils.start("Bootstrapping Helper construct");
	BootHelper boothelper(log2n, radix, logc, scheme, ring, secretKey);
	timeutils.stop("Bootstrapping Helper construct");

    double* mvec = new double[n];
    for(int i = 0; i < n; i++) {
        mvec[i] = (double) rand() / RAND_MAX;
    }
    
	Ciphertext cipher = scheme.encrypt(mvec, n, logp, logQ);

	Ciphertext sortedCipher;
    timeutils.start("encBatcherOddEvenSort");
	encBatcherOddEvenSort(sortedCipher, cipher, parameter, iter, scheme, boothelper);
    timeutils.stop("encBatcherOddEvenSort");

	// run Plain Sorting


	SortingAlgorithm sort(CyclicArray(mvec, 1 << log2n), log2n);
    sort.BatcherOddEvenSort();

	complex<double>* dvec = scheme.decrypt(secretKey, sortedCipher);

	cout <<  "num : Original // PlainSort // EncSort" << endl;
    for(int i = 0; i < n; i++) {
        // cout << i << " : " << mvec[i] << ", " << dvec[i].real() << endl;
		cout << i << " : " << mvec[i] << " // " << sort.ca.get(i) << " // " << dvec[i].real() << endl;
    }
}

void TestBootstrapping::testEncSortWithDecrypt(Parameter parameter, long iter)  {
         // HE parameter //
	long logN = parameter.logN;
	long logQ = parameter.logQ;
	long logp = parameter.logp;
	long logc = parameter.logc;

	// Decomposition related parameter //
	long log2n = parameter.log2n;
	long radix = parameter.radix;

	// Bootstrapping parameter //
	long logq = parameter.logq;
	long logT = parameter.logT;

	long n = 1 << log2n;

	cout << "\n***************************" << endl;
	cout << "Test for EncSort" << endl;
	cout << "logN = " << logN << ", logQ = " << logQ << ", logp = " << logp << ", logc = " << logc << endl;
	cout << "slots = " << n << ", radix = " << radix << ", logq = " << logq << ", logT = " << logT << endl;
	cout << "***************************" << endl;
	cout << endl;

	TimeUtils timeutils;
	timeutils.start("KeyGen");
	Ring ring(logN, logQ);
	SecretKey secretKey(ring);
	Scheme scheme(secretKey, ring);
	scheme.addConjKey(secretKey);
	scheme.addLeftRotKeys(secretKey);
	scheme.addRightRotKeys(secretKey);
	timeutils.stop("KeyGen");
 
	timeutils.start("Bootstrapping Helper construct");
	BootHelper boothelper(log2n, radix, logc, scheme, ring, secretKey);
	timeutils.stop("Bootstrapping Helper construct");

    double* mvec = new double[n];
    for(int i = 0; i < n; i++) {
        mvec[i] = (double) rand() / RAND_MAX;
    }
    
	Ciphertext cipher = scheme.encrypt(mvec, n, logp, logQ);

	Ciphertext sortedCipher;
    timeutils.start("encBatcherOddEvenSort");
	encBatcherOddEvenSortWithDecrypt(sortedCipher, cipher, parameter, iter, scheme, boothelper, secretKey);
    timeutils.stop("encBatcherOddEvenSort");

	// run Plain Sorting


	SortingAlgorithm sort(CyclicArray(mvec, 1 << log2n), log2n);
    sort.BatcherOddEvenSort();

	complex<double>* dvec = scheme.decrypt(secretKey, sortedCipher);

	cout <<  "num : Original // PlainSort // EncSort" << endl;
    for(int i = 0; i < n; i++) {
        // cout << i << " : " << mvec[i] << ", " << dvec[i].real() << endl;
		cout << i << " : " << mvec[i] << " // " << sort.ca.get(i) << " // " << dvec[i].real() << endl;
    }
}

void TestBootstrapping::testSqrtWithBootAndDecrypt(Parameter parameter, long iter) {
	// HE parameter //
	long logN = parameter.logN;
	long logQ = parameter.logQ;
	long logp = parameter.logp;
	long logc = parameter.logc;

	// Decomposition related parameter //
	long log2n = parameter.log2n;
	long radix = parameter.radix;

	// Bootstrapping parameter //
	long logq = parameter.logq;
	long logT = parameter.logT;

	long n = 1 << log2n;

	cout << "\n***************************" << endl;
	cout << "Test for Sqrt with Boot" << endl;
	cout << "logN = " << logN << ", logQ = " << logQ << ", logp = " << logp << ", logc = " << logc << endl;
	cout << "slots = " << n << ", radix = " << radix << ", logq = " << logq << ", logT = " << logT << endl;
	cout << "***************************" << endl;
	cout << endl;

	TimeUtils timeutils;
	timeutils.start("KeyGen");
	Ring ring(logN, logQ);
	SecretKey secretKey(ring);
	Scheme scheme(secretKey, ring);
	scheme.addConjKey(secretKey);
	scheme.addLeftRotKeys(secretKey);
	timeutils.stop("KeyGen");

	timeutils.start("Bootstrapping Helper construct");
	BootHelper boothelper(log2n, radix, logc, scheme, ring, secretKey);
	timeutils.stop("Bootstrapping Helper construct");

	double* mvec = new double[n];
    for(int i = 0; i < n; i++) {
        mvec[i] = 0.00001 + (double) i / n;
    }
    
	Ciphertext cipher = scheme.encrypt(mvec, n, logp, logQ);

    Ciphertext sqrtCipher;
    timeutils.start("Sqrt");
    // fcnEncComputeSqrt(sqrtCipher, cipher, logp, 5, scheme);
	fcnEncSqrtWithBootWithDecrypt(sqrtCipher, cipher, parameter, iter, scheme, boothelper, secretKey);
    timeutils.stop("Sqrt");

    // Print Result and Difference //	
	complex<double>* dvec = scheme.decrypt(secretKey, sqrtCipher);
    for(int i = 0; i < n; i++) {
		double ans = dvec[i].real();
        cout << i << " : " << mvec[i] << " ~simeq~ "
			<< ans * ans << " = " << ans << " ^2" << endl;
    }
}

void TestBootstrapping::testEncCompAndSwapWithBootAndDecrypt(Parameter parameter, long iter) {
         // HE parameter //
	long logN = parameter.logN;
	long logQ = parameter.logQ;
	long logp = parameter.logp;
	long logc = parameter.logc;

	// Decomposition related parameter //
	long log2n = parameter.log2n;
	long radix = parameter.radix;

	// Bootstrapping parameter //
	long logq = parameter.logq;
	long logT = parameter.logT;

	long n = 1 << log2n;

	cout << "\n***************************" << endl;
	cout << "Test for Comparison And Swap with Boot" << endl;
	cout << "logN = " << logN << ", logQ = " << logQ << ", logp = " << logp << ", logc = " << logc << endl;
	cout << "slots = " << n << ", radix = " << radix << ", logq = " << logq << ", logT = " << logT << endl;
	cout << "***************************" << endl;
	cout << endl;

	TimeUtils timeutils;
	timeutils.start("KeyGen");
	Ring ring(logN, logQ);
	SecretKey secretKey(ring);
	Scheme scheme(secretKey, ring);
	scheme.addConjKey(secretKey);
	scheme.addLeftRotKeys(secretKey);
	scheme.addRightRotKeys(secretKey);
	timeutils.stop("KeyGen");

	timeutils.start("Bootstrapping Helper construct");
	BootHelper boothelper(log2n, radix, logc, scheme, ring, secretKey);
	timeutils.stop("Bootstrapping Helper construct");

    double* mvec = new double[n];
    for(int i = 0; i < n; i++) {
        // mvec[n - 1 - i] = 0.00001 + (double) i / n;
		mvec[i] = (double) rand() / RAND_MAX;
    }
    
	Ciphertext cipher = scheme.encrypt(mvec, n, logp, logQ);

	Ciphertext cipher2;

    double* mask = genMaskingComp(n, 1);
    
    timeutils.start("CompAndSwap");
	// fcnEncCompAndSwapWithBootAndDecrypt(dummy, cipher, mask, 1, parameter, iter, scheme, boothelper, secretKey);
	fcnEncCompAndSwapWithBoot(cipher2, cipher, mask, 1, parameter, iter, scheme, boothelper);
    timeutils.stop("CompAndSwap");

	cout << "cipher.logq = " << cipher.logq << endl;
	cout << "cipher2.logq = " << cipher2.logq << endl;

	complex<double>* dvec = scheme.decrypt(secretKey, cipher2);
    // StringUtils::compare(mvec, dmult, n, "mult");
    for(int i = 0; i < n; i++) {
        cout << i << " : " << mvec[i] << ", " << dvec[i].real() << endl;
    }
}