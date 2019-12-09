package NumCalc;

import NumCalc.Properties.*;
import ClassicalStat.*;

//endregion IQrDecomposition

//region Interface ISingularValueDecomposition

/** 
 	Singular Value Decomposition for a rectangular matrix.
 
 
	  For an m-by-n matrix <c>A</c> with <c>m >= n</c>, the singular value decomposition is
   an m-by-n orthogonal matrix <c>U</c>, an n-by-n diagonal matrix <c>S</c>, and
   an n-by-n orthogonal matrix <c>V</c> so that <c>A = U * S * V'</c>.
   The singular values, <c>sigma[df] = S[df,df]</c>, are ordered so that
   <c>sigma[0] >= sigma[1] >= ... >= sigma[n-1]</c>.
   The singular value decompostion always exists, so the constructor will
   never fail. The matrix condition number and the effective numerical
   rank can be computed from this decomposition.
 
*/
public interface ISingularValueDecomposition
{
	/** Returns the condition number <c>max(S) / min(S)</c>.
	*/
	double getCondition();

	/** Returns the Two norm.
	*/
	double getNorm2();

	/** Returns the effective numerical matrix rank.
	 <value>Number of non-negligible singular values.</value>
	*/
	int getRank();

	/** Return the one-dimensional array of singular values.
	*/
	double[] getDiagonal();
}