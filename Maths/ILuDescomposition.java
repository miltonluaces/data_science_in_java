package NumCalc;

import NumCalc.Properties.*;
import ClassicalStat.*;

//endregion ICholeskyDecomposition

//region Interface ILuDecomposition

/** 
   LU decomposition of a rectangular matrix.
 
 
   For an m-by-n matrix <c>A</c> with m >= n, the LU decomposition is an m-by-n
   unit lower triangular matrix <c>L</c>, an n-by-n upper triangular matrix <c>U</c>,
   and a permutation vector <c>piv</c> of length m so that <c>A(piv)=L*U</c>.
   If m &lt; n, then <c>L</c> is m-by-m and <c>U</c> is m-by-n.
   The LU decompostion with pivoting always exists, even if the matrix is
   singular, so the constructor will never fail.  The primary use of the
   LU decomposition is in the solution of square systems of simultaneous
   linear equations. This will fail if <see cref="IsNonSingular"/> returns <see langword="false"/>.
 
*/
public interface ILuDecomposition
{
	/** Returns if the matrix is non-singular.
	*/
	boolean getIsNonSingular();

	/** Returns the determinant of the matrix.
	*/
	double getDeterminant();

	/** Returns the lower triangular factor <c>L</c> with <c>A=LU</c>.
	*/
	IMatrix getLowerTriangularFactor();

	/** Returns the lower triangular factor <c>L</c> with <c>A=LU</c>.
	*/
	IMatrix getUpperTriangularFactor();

	/** Returns the pivot permuation vector.
	*/
	double[] getPivotPermutationVector();

	/** Solves a set of equation systems of type <c>A * X = B</c>.
	 @param rhs Right hand side matrix with as many rows as <c>A</c> and any number of columns.
	 @return Matrix <c>X</c> so that <c>L * U * X = B</c>.
	*/
	IMatrix Solve(IMatrix rhs);
}