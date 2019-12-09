package Clustering;

import ClassicalStat.*;

public class svd
{
	/*************************************************************************
	Singular value decomposition of a rectangular matrix.

	The algorithm calculates the singular value decomposition of a matrix of
	size MxN: A = U * S * V^T

	The algorithm finds the singular values and, optionally, matrices U and V^T.
	The algorithm can find both first min(M,N) columns of matrix U and rows of
	matrix V^T (singular vectors), and matrices U and V^T wholly (of sizes MxM
	and NxN respectively).

	Take into account that the subroutine does not return matrix V but V^T.

	Input parameters:
	    A           -   matrix to be decomposed.
	                    Array whose indexes range within [0..M-1, 0..N-1].
	    M           -   number of rows in matrix A.
	    N           -   number of columns in matrix A.
	    UNeeded     -   0, 1 or 2. See the description of the parameter U.
	    VTNeeded    -   0, 1 or 2. See the description of the parameter VT.
	    AdditionalMemory -
	                    If the parameter:
	                     * equals 0, the algorithm doesn’t use additional
	                       memory (lower requirements, lower performance).
	                     * equals 1, the algorithm uses additional
	                       memory of size min(M,N)*min(M,N) of real numbers.
	                       It often speeds up the algorithm.
	                     * equals 2, the algorithm uses additional
	                       memory of size M*min(M,N) of real numbers.
	                       It allows to get a maximum performance.
	                    The recommended value of the parameter is 2.

	Output parameters:
	    W           -   contains singular values in descending order.
	    U           -   if UNeeded=0, U isn't changed, the left singular vectors
	                    are not calculated.
	                    if Uneeded=1, U contains left singular vectors (first
	                    min(M,N) columns of matrix U). Array whose indexes range
	                    within [0..M-1, 0..Min(M,N)-1].
	                    if UNeeded=2, U contains matrix U wholly. Array whose
	                    indexes range within [0..M-1, 0..M-1].
	    VT          -   if VTNeeded=0, VT isn’t changed, the right singular vectors
	                    are not calculated.
	                    if VTNeeded=1, VT contains right singular vectors (first
	                    min(M,N) rows of matrix V^T). Array whose indexes range
	                    within [0..min(M,N)-1, 0..N-1].
	                    if VTNeeded=2, VT contains matrix V^T wholly. Array whose
	                    indexes range within [0..N-1, 0..N-1].

	  -- ALGLIB --
	     Copyright 2005 by Bochkanov Sergey
	*************************************************************************/
	public static boolean rmatrixsvd(double[][] a, int m, int n, int uneeded, int vtneeded, int additionalmemory, tangible.RefObject<double[]> w, tangible.RefObject<double[][]> u, tangible.RefObject<double[][]> vt)
	{
		boolean result = new boolean();
		double[] tauq = new double[0];
		double[] taup = new double[0];
		double[] tau = new double[0];
		double[] e = new double[0];
		double[] work = new double[0];
		double[][] t2 = new double[0][0];
		boolean isupper = new boolean();
		int minmn = 0;
		int ncu = 0;
		int nrvt = 0;
		int nru = 0;
		int ncvt = 0;
		int i = 0;
		int j = 0;

		a = (double[][])a.clone();
		w.argValue = new double[0];
		u.argValue = new double[0][0];
		vt.argValue = new double[0][0];

		result = true;
		if (m == 0 || n == 0)
		{
			return result;
		}

		//
		// initialize
		//
		minmn = Math.min(m, n);
		w.argValue = new double[minmn + 1];
		ncu = 0;
		nru = 0;
		if (uneeded == 1)
		{
			nru = m;
			ncu = minmn;
			u.argValue = new double[nru - 1 + 1][ncu - 1 + 1];
		}
		if (uneeded == 2)
		{
			nru = m;
			ncu = m;
			u.argValue = new double[nru - 1 + 1][ncu - 1 + 1];
		}
		nrvt = 0;
		ncvt = 0;
		if (vtneeded == 1)
		{
			nrvt = minmn;
			ncvt = n;
			vt.argValue = new double[nrvt - 1 + 1][ncvt - 1 + 1];
		}
		if (vtneeded == 2)
		{
			nrvt = n;
			ncvt = n;
			vt.argValue = new double[nrvt - 1 + 1][ncvt - 1 + 1];
		}

		//
		// M much larger than N
		// Use bidiagonal reduction with QR-decomposition
		//
		if ((double)(m) > (double)(1.6 * n))
		{
			if (uneeded == 0)
			{

				//
				// No left singular vectors to be computed
				//
				tangible.RefObject<Double> tempRef_a = new tangible.RefObject<Double>(a);
				tangible.RefObject<Double> tempRef_tau = new tangible.RefObject<Double>(tau);
				ortfac.rmatrixqr(tempRef_a, m, n, tempRef_tau);
			tau = tempRef_tau.argValue;
			a = tempRef_a.argValue;
				for (i = 0; i <= n - 1; i++)
				{
					for (j = 0; j <= i - 1; j++)
					{
						a[i][j] = 0;
					}
				}
				tangible.RefObject<Double> tempRef_a2 = new tangible.RefObject<Double>(a);
				tangible.RefObject<Double> tempRef_tauq = new tangible.RefObject<Double>(tauq);
				tangible.RefObject<Double> tempRef_taup = new tangible.RefObject<Double>(taup);
				ortfac.rmatrixbd(tempRef_a2, n, n, tempRef_tauq, tempRef_taup);
			taup = tempRef_taup.argValue;
			tauq = tempRef_tauq.argValue;
			a = tempRef_a2.argValue;
				ortfac.rmatrixbdunpackpt(a, n, n, taup, nrvt, vt);
				tangible.RefObject<Boolean> tempRef_isupper = new tangible.RefObject<Boolean>(isupper);
				tangible.RefObject<Double> tempRef_e = new tangible.RefObject<Double>(e);
				ortfac.rmatrixbdunpackdiagonals(a, n, n, tempRef_isupper, w, tempRef_e);
			e = tempRef_e.argValue;
			isupper = tempRef_isupper.argValue;
				tangible.RefObject<double[][]> tempRef_a3 = new tangible.RefObject<double[][]>(a);
				result = bdsvd.rmatrixbdsvd(w, e, n, isupper, false, u, 0, tempRef_a3, 0, vt, ncvt);
			a = tempRef_a3.argValue;
				return result;
			}
			else
			{

				//
				// Left singular vectors (may be full matrix U) to be computed
				//
				tangible.RefObject<Double> tempRef_a4 = new tangible.RefObject<Double>(a);
				tangible.RefObject<Double> tempRef_tau2 = new tangible.RefObject<Double>(tau);
				ortfac.rmatrixqr(tempRef_a4, m, n, tempRef_tau2);
			tau = tempRef_tau2.argValue;
			a = tempRef_a4.argValue;
				ortfac.rmatrixqrunpackq(a, m, n, tau, ncu, u);
				for (i = 0; i <= n - 1; i++)
				{
					for (j = 0; j <= i - 1; j++)
					{
						a[i][j] = 0;
					}
				}
				tangible.RefObject<Double> tempRef_a5 = new tangible.RefObject<Double>(a);
				tangible.RefObject<Double> tempRef_tauq2 = new tangible.RefObject<Double>(tauq);
				tangible.RefObject<Double> tempRef_taup2 = new tangible.RefObject<Double>(taup);
				ortfac.rmatrixbd(tempRef_a5, n, n, tempRef_tauq2, tempRef_taup2);
			taup = tempRef_taup2.argValue;
			tauq = tempRef_tauq2.argValue;
			a = tempRef_a5.argValue;
				ortfac.rmatrixbdunpackpt(a, n, n, taup, nrvt, vt);
				tangible.RefObject<Boolean> tempRef_isupper2 = new tangible.RefObject<Boolean>(isupper);
				tangible.RefObject<Double> tempRef_e2 = new tangible.RefObject<Double>(e);
				ortfac.rmatrixbdunpackdiagonals(a, n, n, tempRef_isupper2, w, tempRef_e2);
			e = tempRef_e2.argValue;
			isupper = tempRef_isupper2.argValue;
				if (additionalmemory < 1)
				{

					//
					// No additional memory can be used
					//
					ortfac.rmatrixbdmultiplybyq(a, n, n, tauq, u, m, n, true, false);
					tangible.RefObject<double[][]> tempRef_a6 = new tangible.RefObject<double[][]>(a);
					result = bdsvd.rmatrixbdsvd(w, e, n, isupper, false, u, m, tempRef_a6, 0, vt, ncvt);
				a = tempRef_a6.argValue;
				}
				else
				{

					//
					// Large U. Transforming intermediate matrix T2
					//
					work = new double[Math.max(m, n) + 1];
					tangible.RefObject<Double> tempRef_t2 = new tangible.RefObject<Double>(t2);
					ortfac.rmatrixbdunpackq(a, n, n, tauq, n, tempRef_t2);
				t2 = tempRef_t2.argValue;
					tangible.RefObject<Double> tempRef_a7 = new tangible.RefObject<Double>(a);
					blas.copymatrix(u.argValue, 0, m - 1, 0, n - 1, tempRef_a7, 0, m - 1, 0, n - 1);
				a = tempRef_a7.argValue;
					tangible.RefObject<Double> tempRef_t22 = new tangible.RefObject<Double>(t2);
					tangible.RefObject<Double> tempRef_work = new tangible.RefObject<Double>(work);
					blas.inplacetranspose(tempRef_t22, 0, n - 1, 0, n - 1, tempRef_work);
				work = tempRef_work.argValue;
				t2 = tempRef_t22.argValue;
					tangible.RefObject<double[][]> tempRef_t23 = new tangible.RefObject<double[][]>(t2);
					result = bdsvd.rmatrixbdsvd(w, e, n, isupper, false, u, 0, tempRef_t23, n, vt, ncvt);
				t2 = tempRef_t23.argValue;
					tangible.RefObject<Double> tempRef_work2 = new tangible.RefObject<Double>(work);
					blas.matrixmatrixmultiply(a, 0, m - 1, 0, n - 1, false, t2, 0, n - 1, 0, n - 1, true, 1.0, u, 0, m - 1, 0, n - 1, 0.0, tempRef_work2);
				work = tempRef_work2.argValue;
				}
				return result;
			}
		}

		//
		// N much larger than M
		// Use bidiagonal reduction with LQ-decomposition
		//
		if ((double)(n) > (double)(1.6 * m))
		{
			if (vtneeded == 0)
			{

				//
				// No right singular vectors to be computed
				//
				tangible.RefObject<Double> tempRef_a8 = new tangible.RefObject<Double>(a);
				tangible.RefObject<Double> tempRef_tau3 = new tangible.RefObject<Double>(tau);
				ortfac.rmatrixlq(tempRef_a8, m, n, tempRef_tau3);
			tau = tempRef_tau3.argValue;
			a = tempRef_a8.argValue;
				for (i = 0; i <= m - 1; i++)
				{
					for (j = i + 1; j <= m - 1; j++)
					{
						a[i][j] = 0;
					}
				}
				tangible.RefObject<Double> tempRef_a9 = new tangible.RefObject<Double>(a);
				tangible.RefObject<Double> tempRef_tauq3 = new tangible.RefObject<Double>(tauq);
				tangible.RefObject<Double> tempRef_taup3 = new tangible.RefObject<Double>(taup);
				ortfac.rmatrixbd(tempRef_a9, m, m, tempRef_tauq3, tempRef_taup3);
			taup = tempRef_taup3.argValue;
			tauq = tempRef_tauq3.argValue;
			a = tempRef_a9.argValue;
				ortfac.rmatrixbdunpackq(a, m, m, tauq, ncu, u);
				tangible.RefObject<Boolean> tempRef_isupper3 = new tangible.RefObject<Boolean>(isupper);
				tangible.RefObject<Double> tempRef_e3 = new tangible.RefObject<Double>(e);
				ortfac.rmatrixbdunpackdiagonals(a, m, m, tempRef_isupper3, w, tempRef_e3);
			e = tempRef_e3.argValue;
			isupper = tempRef_isupper3.argValue;
				work = new double[m + 1];
				tangible.RefObject<Double> tempRef_work3 = new tangible.RefObject<Double>(work);
				blas.inplacetranspose(u, 0, nru - 1, 0, ncu - 1, tempRef_work3);
			work = tempRef_work3.argValue;
				tangible.RefObject<double[][]> tempRef_a10 = new tangible.RefObject<double[][]>(a);
				result = bdsvd.rmatrixbdsvd(w, e, m, isupper, false, tempRef_a10, 0, u, nru, vt, 0);
			a = tempRef_a10.argValue;
				tangible.RefObject<Double> tempRef_work4 = new tangible.RefObject<Double>(work);
				blas.inplacetranspose(u, 0, nru - 1, 0, ncu - 1, tempRef_work4);
			work = tempRef_work4.argValue;
				return result;
			}
			else
			{

				//
				// Right singular vectors (may be full matrix VT) to be computed
				//
				tangible.RefObject<Double> tempRef_a11 = new tangible.RefObject<Double>(a);
				tangible.RefObject<Double> tempRef_tau4 = new tangible.RefObject<Double>(tau);
				ortfac.rmatrixlq(tempRef_a11, m, n, tempRef_tau4);
			tau = tempRef_tau4.argValue;
			a = tempRef_a11.argValue;
				ortfac.rmatrixlqunpackq(a, m, n, tau, nrvt, vt);
				for (i = 0; i <= m - 1; i++)
				{
					for (j = i + 1; j <= m - 1; j++)
					{
						a[i][j] = 0;
					}
				}
				tangible.RefObject<Double> tempRef_a12 = new tangible.RefObject<Double>(a);
				tangible.RefObject<Double> tempRef_tauq4 = new tangible.RefObject<Double>(tauq);
				tangible.RefObject<Double> tempRef_taup4 = new tangible.RefObject<Double>(taup);
				ortfac.rmatrixbd(tempRef_a12, m, m, tempRef_tauq4, tempRef_taup4);
			taup = tempRef_taup4.argValue;
			tauq = tempRef_tauq4.argValue;
			a = tempRef_a12.argValue;
				ortfac.rmatrixbdunpackq(a, m, m, tauq, ncu, u);
				tangible.RefObject<Boolean> tempRef_isupper4 = new tangible.RefObject<Boolean>(isupper);
				tangible.RefObject<Double> tempRef_e4 = new tangible.RefObject<Double>(e);
				ortfac.rmatrixbdunpackdiagonals(a, m, m, tempRef_isupper4, w, tempRef_e4);
			e = tempRef_e4.argValue;
			isupper = tempRef_isupper4.argValue;
				work = new double[Math.max(m, n) + 1];
				tangible.RefObject<Double> tempRef_work5 = new tangible.RefObject<Double>(work);
				blas.inplacetranspose(u, 0, nru - 1, 0, ncu - 1, tempRef_work5);
			work = tempRef_work5.argValue;
				if (additionalmemory < 1)
				{

					//
					// No additional memory available
					//
					ortfac.rmatrixbdmultiplybyp(a, m, m, taup, vt, m, n, false, true);
					tangible.RefObject<double[][]> tempRef_a13 = new tangible.RefObject<double[][]>(a);
					result = bdsvd.rmatrixbdsvd(w, e, m, isupper, false, tempRef_a13, 0, u, nru, vt, n);
				a = tempRef_a13.argValue;
				}
				else
				{

					//
					// Large VT. Transforming intermediate matrix T2
					//
					tangible.RefObject<Double> tempRef_t24 = new tangible.RefObject<Double>(t2);
					ortfac.rmatrixbdunpackpt(a, m, m, taup, m, tempRef_t24);
				t2 = tempRef_t24.argValue;
					tangible.RefObject<double[][]> tempRef_a14 = new tangible.RefObject<double[][]>(a);
					tangible.RefObject<double[][]> tempRef_t25 = new tangible.RefObject<double[][]>(t2);
					result = bdsvd.rmatrixbdsvd(w, e, m, isupper, false, tempRef_a14, 0, u, nru, tempRef_t25, m);
				t2 = tempRef_t25.argValue;
				a = tempRef_a14.argValue;
					tangible.RefObject<Double> tempRef_a15 = new tangible.RefObject<Double>(a);
					blas.copymatrix(vt.argValue, 0, m - 1, 0, n - 1, tempRef_a15, 0, m - 1, 0, n - 1);
				a = tempRef_a15.argValue;
					tangible.RefObject<Double> tempRef_work6 = new tangible.RefObject<Double>(work);
					blas.matrixmatrixmultiply(t2, 0, m - 1, 0, m - 1, false, a, 0, m - 1, 0, n - 1, false, 1.0, vt, 0, m - 1, 0, n - 1, 0.0, tempRef_work6);
				work = tempRef_work6.argValue;
				}
				tangible.RefObject<Double> tempRef_work7 = new tangible.RefObject<Double>(work);
				blas.inplacetranspose(u, 0, nru - 1, 0, ncu - 1, tempRef_work7);
			work = tempRef_work7.argValue;
				return result;
			}
		}

		//
		// M<=N
		// We can use inplace transposition of U to get rid of columnwise operations
		//
		if (m <= n)
		{
			tangible.RefObject<Double> tempRef_a16 = new tangible.RefObject<Double>(a);
			tangible.RefObject<Double> tempRef_tauq5 = new tangible.RefObject<Double>(tauq);
			tangible.RefObject<Double> tempRef_taup5 = new tangible.RefObject<Double>(taup);
			ortfac.rmatrixbd(tempRef_a16, m, n, tempRef_tauq5, tempRef_taup5);
		taup = tempRef_taup5.argValue;
		tauq = tempRef_tauq5.argValue;
		a = tempRef_a16.argValue;
			ortfac.rmatrixbdunpackq(a, m, n, tauq, ncu, u);
			ortfac.rmatrixbdunpackpt(a, m, n, taup, nrvt, vt);
			tangible.RefObject<Boolean> tempRef_isupper5 = new tangible.RefObject<Boolean>(isupper);
			tangible.RefObject<Double> tempRef_e5 = new tangible.RefObject<Double>(e);
			ortfac.rmatrixbdunpackdiagonals(a, m, n, tempRef_isupper5, w, tempRef_e5);
		e = tempRef_e5.argValue;
		isupper = tempRef_isupper5.argValue;
			work = new double[m + 1];
			tangible.RefObject<Double> tempRef_work8 = new tangible.RefObject<Double>(work);
			blas.inplacetranspose(u, 0, nru - 1, 0, ncu - 1, tempRef_work8);
		work = tempRef_work8.argValue;
			tangible.RefObject<double[][]> tempRef_a17 = new tangible.RefObject<double[][]>(a);
			result = bdsvd.rmatrixbdsvd(w, e, minmn, isupper, false, tempRef_a17, 0, u, nru, vt, ncvt);
		a = tempRef_a17.argValue;
			tangible.RefObject<Double> tempRef_work9 = new tangible.RefObject<Double>(work);
			blas.inplacetranspose(u, 0, nru - 1, 0, ncu - 1, tempRef_work9);
		work = tempRef_work9.argValue;
			return result;
		}

		//
		// Simple bidiagonal reduction
		//
		tangible.RefObject<Double> tempRef_a18 = new tangible.RefObject<Double>(a);
		tangible.RefObject<Double> tempRef_tauq6 = new tangible.RefObject<Double>(tauq);
		tangible.RefObject<Double> tempRef_taup6 = new tangible.RefObject<Double>(taup);
		ortfac.rmatrixbd(tempRef_a18, m, n, tempRef_tauq6, tempRef_taup6);
	taup = tempRef_taup6.argValue;
	tauq = tempRef_tauq6.argValue;
	a = tempRef_a18.argValue;
		ortfac.rmatrixbdunpackq(a, m, n, tauq, ncu, u);
		ortfac.rmatrixbdunpackpt(a, m, n, taup, nrvt, vt);
		tangible.RefObject<Boolean> tempRef_isupper6 = new tangible.RefObject<Boolean>(isupper);
		tangible.RefObject<Double> tempRef_e6 = new tangible.RefObject<Double>(e);
		ortfac.rmatrixbdunpackdiagonals(a, m, n, tempRef_isupper6, w, tempRef_e6);
	e = tempRef_e6.argValue;
	isupper = tempRef_isupper6.argValue;
		if (additionalmemory < 2 || uneeded == 0)
		{

			//
			// We cant use additional memory or there is no need in such operations
			//
			tangible.RefObject<double[][]> tempRef_a19 = new tangible.RefObject<double[][]>(a);
			result = bdsvd.rmatrixbdsvd(w, e, minmn, isupper, false, u, nru, tempRef_a19, 0, vt, ncvt);
		a = tempRef_a19.argValue;
		}
		else
		{

			//
			// We can use additional memory
			//
			t2 = new double[minmn - 1 + 1][m - 1 + 1];
			tangible.RefObject<Double> tempRef_t26 = new tangible.RefObject<Double>(t2);
			blas.copyandtranspose(u.argValue, 0, m - 1, 0, minmn - 1, tempRef_t26, 0, minmn - 1, 0, m - 1);
		t2 = tempRef_t26.argValue;
			tangible.RefObject<double[][]> tempRef_t27 = new tangible.RefObject<double[][]>(t2);
			result = bdsvd.rmatrixbdsvd(w, e, minmn, isupper, false, u, 0, tempRef_t27, m, vt, ncvt);
		t2 = tempRef_t27.argValue;
			blas.copyandtranspose(t2, 0, minmn - 1, 0, m - 1, u, 0, m - 1, 0, minmn - 1);
		}
		return result;
	}
}