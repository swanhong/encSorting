#ifndef NEW_BOOTSTRAPPING_H_
#define NEW_BOOTSTRAPPING_H_

#include "common.h"

class BootHelper
{

public:

	// ring class for encoding
	Ring& ring;

	// scheme class for homomorphic computation
	Scheme& scheme;

	// parameters
	long slots;
	long radix;
	long logSlots;
	long logrSlots;
	long log2r;
	long logp;

	// decomposed linear transformation in bootstrapping
	MatrixXc* V0Decomp;
	MatrixXc* invV0Decomp;

	// pre-encoded polys
	ZZ*** encodeV0poly;
	ZZ*** encodeinvV0poly;
	
	BootHelper(long logSlots, long radix, long logp, Scheme& scheme, Ring& ring, SecretKey& secretKey);

	~BootHelper();

	void buildPartialMatrix(MatrixXc& A, long start_row, long start_col, long size);

	void buildV0DecompMatrix(MatrixXc* DFTDecomp, long n);

	void encodeDiagonal(ZZ* poly, MatrixXc A, long idx);

	void V0(Ciphertext& cipher);

	void invV0(Ciphertext& cipher);

	void coeffToSlot(Ciphertext& part1, Ciphertext& part2, Ciphertext& cipher);

	void slotToCoeff(Ciphertext& cipher, Ciphertext& part1, Ciphertext& part2);

	void evalExpAndEqual(Ciphertext& part1, Ciphertext& part2, long logT, long logI, long logq);

	void bootstrapping(Ciphertext& cipher, long logq, long logQ, long logT);

};


#endif // ! NEW_BOOTSTRAPPING_H_