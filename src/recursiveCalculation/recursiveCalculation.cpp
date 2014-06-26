/*
 * recursiveCalculation.cpp
 *
 *  Created on: Jun 24, 2014
 *      Author: pxiang
 */

#include "recursiveCalculation.h"

/**
 * clean up binary files that are used to save the matrices
 */
void deleteMatrixFiles(std::string files) {
	std::string command = "rm " + files;
	system(command.c_str());
}

/**
 * solve the linear equation A*X = B
 */
void solveDenseLinearEqs(CDMatrix& A, CDMatrix& B, CDMatrix& X) {
	X = A.colPivHouseholderQr().solve(B);
}

/**
 * set up the precondition for the recursive calculations: fromRightToCenter
 * and fromLeftToCenter
 */
void setUpRecursion(LatticeShape& lattice, Parameters& parameters,
		            Basis& initialSites, int maxDistance,
		            int& KLeftStart, int& KLeftStop, int& KCritial,
		            int& KRightStop, int& KRightStart) {
	switch (lattice.getDim()) {
	case 1:
		int Kmin = 1;
		int Kmax = 2*lattice.getXmax() - 1;
		int Kc = initialSites[0] + initialSites[1];

		KLeftStart = 1; // always start from V_1 from the left
		// find out the position where the left and right recursions must stop
		for (int K=1; K<=Kmax; K+=maxDistance) {
			if (Kc>=K && Kc<=K+maxDistance-1) {
				KCritial = K;
				KLeftStop = KCritial-maxDistance;
				KRightStop = KCritial+maxDistance;
			}
		}

		/**
		 * To find out KRightStart, we start from
		 *  V_{KRightStop}=[ v_{KRightStop}, v_{KRightStop+1},
		 *                  ..., v_{KRightStop+maxDistance-1} ].
		 * We can let K = KRightStop+maxDistance-1 and increase K by maxDistance
		 * each time until K > Kmax, then the value of KRightStart is given by
		 * current_value_of_K - maxDistance - (maxDistance-1)
		 *
		 */
		int K= KRightStop + maxDistance -1;
		while (K <= Kmax) {
			K += maxDistance;
		}
		KRightStart = K - maxDistance - (maxDistance-1);

		// calculate the indexMatrix and set up the interaction matrix
		generateIndexMatrix(lattice);
		setInteractions(lattice, parameters);
		break;
	case 2:
		break;
	case 3:
		break;
	}
}




/**
 * recursive calculation from right boundary to the center
 *
 * must call setUpRecursion before calling this function
 */
void fromRightToCenter(int KRightStart, int KRightStop, int maxDistance,
		dcomplex z, CDMatrix& AKRightStop, bool saveAMatrices) {
	/**
	 * The recursive relation is given by:
	 * W_{K}*V_{K} = alpha_{K}*V_{K-maxDistance} + beta_{K}*V_{K+maxDistance}
	 * where V_{K} = [v_{K}, v_{K+1}, ..., v_{K+maxDistance-1}] and v_{K} is
	 * composed of the Green's functions whose parameters adds up to K (for
	 * example in 1D, v_{5} is [G(0,5), G(1,4), G(2, 3)] )
	 *
	 * At the right-most end, the index K = KRightStart. But KRightStart is not
	 * necessarily equal to Kmax (the largest value that K can take).
	 * What we can say is that Kmax is one of {KRightStop, KRightStop+1, ...,
	 * KRightStop+maxDistance-1}. Similarly, KRightStop-1 is not necessarily
	 * equal to the summation of initial indexes. For example, in 1D, we
	 * want to calculate G(n, m, n', m'), Kc = n' + m', but KRightStop may be
	 * != Kc + 1. We can only say that V_{KRightStop-1} contains v_{Kc} (that
	 * is why the recursive calculation from right must stop at KRightStop)
	 *
	 * For the convenience of calculation, we may let the Kmax to be
	 * (N*maxDistance) so that all v_{K} from v_1 to v_Kmax can be divided into
	 * V_{1} = [v_1, v_2, ..., v_{maxDistance}],
	 * V_{maxDistance+1} = [v_{maxDistance+1}, v_{maxDistance+2}, ..., v_{2*maxDistance}],
	 * V_{2*maxDistance+1} = [..., V_{3*maxDistance}]
	 *  .
	 *  .
	 *  .
	 * V_{(N-1)*maxDistance+1} = [..., V_{N*maxDistance}]
	 *
	 *
	 * Assuming V_{K+maxDistance} = 0 in the recursive relation for
	 *  K=KRightStart, we obtain:
	 *  W_{KRightStart}*V_{KRightStart} = alpha_{KRightStart}*V_{KRightStart-maxDistance}
	 * Based on the definition: V_{K} = A_{K} * V_{K-maxDistance}, we obtain:
	 * W_{KRightStart} * A_{KRightStart} = alpha_{KRightStart}
	 *
	 * This equation can be solved to give A_{KRightStart}
	 */
	CDMatrix alphaStart;
	formMatrixAlpha(KRightStart,maxDistance, alphaStart);
	/**
	 * W_{K}*V_{K} = alpha_{K}*V_{K-maxDistance} + beta_{K}*V_{K+maxDistance}
	 * let's simplify the notation:
	 * W_{K}*V_{K} = alpha_{K}*V_{K-} + beta_{K}*V_{K+}
	 */
	CDMatrix WKPlus; //initially equal to W_{KRightStart}
	formMatrixW(KRightStart, maxDistance, z, WKPlus);

	CDMatrix AKPlus; //initially equal to A_{KRightStart}
	solveDenseLinearEqs(WKPlus, alphaStart, AKPlus);

	// release memory
	alphaStart.resize(0,0);
	WKPlus.resize(0,0);

	// save the A matrix into a binary file
	std::string filename;
	if (saveAMatrices==true) {
		filename="A"+ itos(KRightStart) + ".bin";
		saveMatrix(filename, AKPlus);
	}

	/**
	 * Now knowing the A_{KRightStart} or AKPlus, we can recursively calculate
	 * AK, AKMinus, ... until AKRightStop
	 *
	 * Start with the recursion relation
	 * W_{K}*V_{K} = alpha_{K}*V_{K-} + beta_{K}*V_{K+}
	 *
	 * Once we know A_{K+}, we can express V_{K+} in terms of V_{K}, that is,
	 * V_{K+} = A_{K+}*V_{K}. Substituting this into the above recursion relation,
	 * we obtain:
	 *  W_{K}*V_{K} = alpha_{K}*V_{K-} + beta_{K}*A_{K+}*V_{K}.
	 *
	 * Because of V_{K} = A_{K} V_{K-}, the above equation can be rewritten as
	 *  W_{K}*A_{K}*V_{K-} = alpha_{K}*V_{K-} + beta_{K}*A_{K+}*A_{K}*V_{K} .
	 * Therefore,
	 * 	W_{K}*A_{K} = alpha_{K} + beta_{K}*A_{K+}*A_{K}
	 *
	 * 	==> ( W_{K} - beta_{K}*A_{K+} ) * A_{K} = alpha_{K}
	 *
	 * 	Then A_{K} can be obtained by solving the above linear equations
	 */
//	CDMatrix tmp1;
//
	for (int K=KRightStart-maxDistance; K>=KRightStop; K-=maxDistance) {
		CDMatrix BetaK;
		formMatrixBeta(K, maxDistance, BetaK);
		CDMatrix WK;
		formMatrixW(K, maxDistance,z, WK);
		CDMatrix LeftSide = WK - BetaK*AKPlus;

		//release the memory of WK and Beta
		BetaK.resize(0,0);
		WK.resize(0,0);

		CDMatrix AlphaK;
		formMatrixAlpha(K, maxDistance, AlphaK);

		// solve for AK and assign the value to AKPlus for next iteration
		solveDenseLinearEqs(LeftSide, AlphaK, AKPlus);

		//release the memory of LeftSide and AlphaK
		LeftSide.resize(0,0);
		AlphaK.resize(0,0);

		// save the AK matrix into a binary file
		if (saveAMatrices==true) {
			filename="A"+ itos(K) + ".bin";
			saveMatrix(filename, AKPlus);
		}
	}

	AKRightStop = AKPlus;
	AKPlus.resize(0,0);
}





/**
 * recursive calculation from the left boundary to the center
 *
 * must call setUpRecursion before calling this function
 */
void fromLeftToCenter(int KLeftStart, int KLeftStop, int maxDistance,
		dcomplex z, CDMatrix& ATildeKLeftStop, bool saveAMatrices) {
	/**
	 * The recursive relation is given by:
	 * W_{K}*V_{K} = alpha_{K}*V_{K-maxDistance} + beta_{K}*V_{K+maxDistance}
	 * where V_{K} = [v_{K}, v_{K+1}, ..., v_{K+maxDistance-1}] and v_{K} is
	 * composed of the Green's functions whose parameters adds up to K (for
	 * example in 1D, v_{5} is [G(0,5), G(1,4), G(2, 3)] )
	 *
	 * At the left-most end, the index K = KLeftStart = Kmin = 1;
	 *
	 *
	 * V_{1} = [v_1, v_2, ..., v_{maxDistance}],
	 * V_{maxDistance+1} = [v_{maxDistance+1}, v_{maxDistance+2}, ..., v_{2*maxDistance}],
	 * V_{2*maxDistance+1} = [..., V_{3*maxDistance}]
	 *  .
	 *  .
	 *  .
	 * V_{(N-1)*maxDistance+1} = [..., V_{N*maxDistance}]
	 *
	 *
	 * Assuming V_{K-maxDistance} = 0 in the recursive relation for
	 * K=KLeftStart, we obtain:
	 * W_{KLeftStart}*V_{KLeftStart} = beta_{KLeftStart}*V_{KLeftStart+maxDistance}
	 * Based on the definition: V_{K} = ATilde_{K} * V_{K+maxDistance}, we obtain:
	 * W_{KLeftStart} * ATilde_{KLeftStart} = beta_{KLeftStart}
	 *
	 * This equation can be solved to give ATilde_{KLeftStart}
	 */
	CDMatrix betaStart;
	formMatrixBeta(KLeftStart,maxDistance, betaStart);
	/**
	 * W_{K}*V_{K} = alpha_{K}*V_{K-maxDistance} + beta_{K}*V_{K+maxDistance}
	 * let's simplify the notation:
	 * W_{K}*V_{K} = alpha_{K}*V_{K-} + beta_{K}*V_{K+}
	 */
	CDMatrix WKMinus; //initially equal to W_{KLeftStart}
	formMatrixW(KLeftStart, maxDistance, z, WKMinus);

	CDMatrix ATildeKMinus; //initially equal to ATilde_{KLeftStart}
	solveDenseLinearEqs(WKMinus, betaStart, ATildeKMinus);

	// release memory
	betaStart.resize(0,0);
	WKMinus.resize(0,0);

	// save the ATilde matrix into a binary file
	std::string filename;
	if (saveAMatrices==true) {
		filename="ATilde"+ itos(KLeftStart) + ".bin";
		saveMatrix(filename, ATildeKMinus);
	}

	/**
	 * Now knowing the A_{KLeftStart} or ATildeKMinus, we can recursively calculate
	 * ATilde, ... until ATildeKLeftStop
	 *
	 * Start with the recursion relation
	 * W_{K}*V_{K} = alpha_{K}*V_{K-} + beta_{K}*V_{K+}
	 *
	 * Once we know ATilde_{K-}, we can express V_{K-} in terms of V_{K},
	 * that is, V_{K-} = ATilde_{K-}*V_{K}. Substituting this into the above
	 * recursion relation, we obtain:
	 *  W_{K}*V_{K} = alpha_{K}*ATilde_{K-}*V_{K} + beta_{K}*V_{K+}.
	 *
	 * Because of V_{K} = ATilde_{K} V_{K+}, the above equation can be rewritten as
	 *  W_{K}*ATilde_{K} V_{K+} = alpha_{K}*ATilde_{K-}*ATilde_{K} V_{K+}
	 *                            + beta_{K}*A_{K+} .
	 * Therefore,
	 * 	W_{K}*ATilde_{K} = alpha_{K}*ATilde_{K-}*ATilde_{K} + beta_{K}
	 *
	 * 	==> ( W_{K} - alpha_{K}*ATilde_{K-} ) * ATilde_{K} = beta_{K}
	 *
	 * 	Then ATilde_{K} can be obtained by solving the above linear equations
	 */
//
	for (int K=KLeftStart+maxDistance; K<=KLeftStop; K += maxDistance) {
		CDMatrix AlphaK;
		formMatrixAlpha(K, maxDistance, AlphaK);
		CDMatrix WK;
		formMatrixW(K, maxDistance,z, WK);
		CDMatrix LeftSide = WK - AlphaK*ATildeKMinus;

		//release the memory of WK and AlphaK
		AlphaK.resize(0,0);
		WK.resize(0,0);

		CDMatrix BetaK;
		formMatrixBeta(K, maxDistance, BetaK);

		// solve for ATildeK and assign the value to ATildeKMinus for next iteration
		solveDenseLinearEqs(LeftSide, BetaK, ATildeKMinus);

		//release the memory of LeftSide and BetaK
		LeftSide.resize(0,0);
		BetaK.resize(0,0);

		// save the AK matrix into a binary file
		if (saveAMatrices==true) {
			filename="ATilde"+ itos(K) + ".bin";
			saveMatrix(filename, ATildeKMinus);
		}
	}

	ATildeKLeftStop = ATildeKMinus;
	ATildeKMinus.resize(0,0);
}


//
//
//
////CDMatrix solveVnc(int ni1, int ni2, complex_double z, Parameters& pars) {
////	AlphaBeta ab(pars);
//CDMatrix solveVnc(int ni1, int ni2, complex_mkl z, AlphaBeta& ab, bool saved) {
//	//	AlphaBeta ab(pars);
//	// to speed up, we can make ab static --- but it doesn't work
//	//static AlphaBeta ab(pars);
//	int Kc = ni1 + ni2;
//	CDMatrix alphanc;
//	CDMatrix betanc;
//	ab.FillAlphaBetaMatrix(Kc,z,alphanc, betanc);
//	CDMatrix AncPlusOne = fromRightToCenter(Kc,z,ab, saved);
//	CDMatrix AncMinusOneTilde = fromLeftToCenter(Kc,z,ab, saved);
//	CDMatrix mat1 = alphanc*AncMinusOneTilde;
//	mat1.AddInPlace(betanc*AncPlusOne);
//	changeElements(mat1);
//
//	// calculate the dimension of the const vector C
//	int nsum = Kc;
//	int nth = 0;
//	int m;
//	for (int i=0; i<= nsum/2; ++i) {
//		int j = nsum - i; // i+j = nsum
//		if (j>i && j<=ab.GetNmax()) {
//			if (i==ni1 && j==ni2) {
//				m = nth;
//			}
//			nth++;
//		}
//	}
//
//
//	CDMatrix constVector(nth, 1, true);
//	complex_mkl denominator =  ab.GetFactor(ni1, ni2, z);
//	//constVector(m,0) = complex_double(1.0, 0.0)/denominator;
//	double square = denominator.real*denominator.real + denominator.imag*denominator.imag;
//	constVector(m,0).real = denominator.real/square;
//	constVector(m,0).imag = -denominator.imag/square;
//
//	solveLinearEqs(mat1, constVector); // now constVector contains the solution
//	return constVector;
//}
//
//

//
//
//
//void generateDensityOfState(int ni1, int ni2, Parameters& pars, const std::vector<complex_mkl>& zList, std::vector<double>& rhoList) {
//	complex_mkl z;
//	CDMatrix Vnc;
//	AlphaBeta ab(pars);
//
//	int nth;
//	rhoList.empty();
//	for (int i=0; i<zList.size(); ++i) {
//		z = zList[i];
//		Vnc = solveVnc(ni1,ni2,z,ab);
//		nth = getIndex(pars.nmax, ni1+ni2, ni1, ni2);
//		rhoList.push_back(-Vnc(nth,0).imag/M_PI);
//	}
//}
//
//

//
//// calculate all the matrix elements of the Green function < i, j | G | m, n> where m and n are fixed
//void calculateAllGF(int ni1, int ni2, complex_mkl z, AlphaBeta& ab) {
//	int nc = ni1 + ni2;
//	std::string fileReadFrom, fileWriteTo;
//	CDMatrix V_K = solveVnc(ni1, ni2, z, ab, true);
//	fileWriteTo="V_"+ itos(nc) + ".bin";
//	CDMatrixToBytes(V_K,fileWriteTo);
//	CDMatrix A;
//	CDMatrix ATilde;
//	CDMatrix V_nc_save = V_K;
//
//	// calculate all V_K and save them into binary files
//	int nmax = ab.GetNmax();
//	// calculate V_K where nc+1 < K <= nmax+nmax-1;
//	for (int K=nc+1; K<=nmax+nmax-1; ++K) {
//		fileReadFrom="A"+ itos(K) + ".bin";
//		bytesToCDMatrix(A,fileReadFrom);
//		V_K = A*V_K;
//		fileWriteTo="V_"+ itos(K) + ".bin";
//		CDMatrixToBytes(V_K,fileWriteTo);
//	}
//
//	// calculate V_K where 1 <= K <= nc-1;
//	V_K = V_nc_save; //restore the value of V_nc
//	for (int K=nc-1; K>=1; --K) {
//		fileReadFrom="ATilde"+ itos(K) + ".bin";
//		bytesToCDMatrix(ATilde,fileReadFrom);
//		V_K = ATilde*V_K;
//		fileWriteTo="V_"+ itos(K) + ".bin";
//		CDMatrixToBytes(V_K,fileWriteTo);
//	}
//}
//
//
//// extract the matrix element G(n, m, ni1, ni2) from files stored in disk
//// it requires the index matrix which tells the order of G(n, m, ni1, ni2) in V_{n+m}
//complex_mkl extractMatrixElement(int n, int m, int ni1, int ni2, IMatrix& indexMatrix) {
//	int nsum = n + m;
//	std::string filename;
//	filename = "V_" + itos(nsum) + ".bin";
//	CDMatrix V_nsum;
//	bytesToCDMatrix(V_nsum, filename);
//	int nth = indexMatrix(n, m);
//	return V_nsum(nth, 0);
//}
