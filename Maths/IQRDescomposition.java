package NumCalc;

import NumCalc.Properties.*;
import ClassicalStat.*;

//endregion ILuDecomposition

//region Interface IQrDecomposition

/** 
	  QR decomposition for a rectangular matrix.
 
 
   For an m-by-n matrix <c>A</c> with <c>m &gt;= n</c>, the QR decomposition is an m-by-n
   orthogonal matrix <c>Q</c> and an n-by-n upper triangular 
   matrix <c>R</c> so that <c>A = Q * R</c>.
   The QR decompostion always exists, even if the matrix does not have
   full rank, so the constructor will never fail.  The primary use of the
   QR decomposition is in the least squares solution of nonsquare systems
   of simultaneous linear equations.
   This will fail if <see cref="IsFullRank"/> returns <see langword="false"/>.
 
*/
public interface IQrDecomposition
{
	/** Shows if the matrix <c>A</c> is of full rank.
	 <value>The value is <see langword="true"/> if <c>R</c>, and hence <c>A</c>, has full rank.</value>
	*/
	boolean getIsFullRank();

	/** Returns the upper triangular factor <c>R</c>.
	*/
	IMatrix getUpperTriangularFactor();

	/** Returns the orthogonal factor <c>Q</c>.
	*/
	IMatrix getOrthogonalFactor();

	/** Least squares solution of <c>A * X = B</c>
	 @param rhs Right-hand-side matrix with as many rows as <c>A</c> and any number of columns.
	 @return A matrix that minimized the two norm of <c>Q * R * X - B</c>.
	 @exception T:System.ArgumentException Matrix row dimensions must be the same.
	 @exception T:System.InvalidOperationException Matrix is rank deficient.
	*/
	IMatrix Solve(IMatrix rhs);
}