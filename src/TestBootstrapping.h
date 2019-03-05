// #include "../bootsrc/new_bootstrapping.h"
#include "Parameter.h"
#include "BootContext.h"
#include "stdlib.h"


class TestBootstrapping {
public:
    static void bootstrapping_test(Parameter parameter);
	
	static void testSqrtWithBoot(Parameter parameter, long iter);
	
	static void testMaxMinWithBoot(Parameter parameter, long iter);

	static void testEncCompAndSwapWithBoot(Parameter parameter, long iter);

	static void testEncSort(Parameter parameter, long iter);

	static void testEncSortWithDecrypt(Parameter parameter, long iter);
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
	Ciphertext cipher = scheme.encrypt(mvec, n, logp, logq);
	mvec = scheme.decrypt(secretKey, cipher);

	complex<double>* dvec;
	
	// Bootstrapping //
	timeutils.start("Improved bootstrapping");
	boothelper.bootstrapping(cipher, logq, logQ, logT);
	timeutils.stop("Improved bootstrapping");

	cout << "* Befor logQ = " << logq << endl;
	cout << "* After logQ = " << cipher.logq << endl;

	// Print Result and Difference //
	dvec = scheme.decrypt(secretKey, cipher);
	cout << "log2(avg of error) = " << diff(mvec, dvec, n) << endl;

	if(mvec != NULL) delete[] mvec;
	if(dvec != NULL) delete[] dvec;
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
        mvec[n - 1 - i] = 0.00001 + (double) i / n;
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