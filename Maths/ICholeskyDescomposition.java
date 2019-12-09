package NumCalc;

import NumCalc.Properties.*;
import ClassicalStat.*;

//endregion IMatrix

//region Interface ICholeskyDecomposition

/** 
		Cholesky Decomposition of a symmetric, positive definite matrix.
	
 
		For a symmetric, positive definite matrix <c>A</c>, the Cholesky decomposition is a
		lower triangular matrix <c>L</c> so that <c>A = L * L'</c>.
		If the matrix is not symmetric or positive definite, the constructor returns a partial 
		decomposition and sets two internal variables that can be queried using the
		<see cref="IsSymmetric"/> and <see cref="IsPositiveDefinite"/> properties.
	
*/
public interface ICholeskyDecomposition
{
	/** Solves a set of equation systems of type <c>A * X = B</c>.
	 @param rhs Right hand side matrix with as many rows as <c>A</c> and any number of columns.
	 @return Matrix <c>X</c> so that <c>L * L' * X = B</c>.
	 @exception T:System.ArgumentException Matrix dimensions do not match.
	 @exception T:System.InvalidOperationException Matrix is not symmetrix and positive definite.
	*/
	IMatrix Solve(IMatrix rhs);

	/** Returns <see langword="true"/> if the matrix is positive definite.
	*/
	boolean getIsPositiveDefinite();

	/** Returns <see langword="true"/> if the matrix is symmetric.
	*/
	boolean getIsSymmetric();

	/** Returns the left triangular factor <c>L</c> so that <c>A = L * L'</c>.
	*/
	IMatrix getLeftTriangularFactor();
}