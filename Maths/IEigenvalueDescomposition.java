package NumCalc;


import ClassicalStat.*;

//endregion ISingularValueDecomposition

//region Interface IEigenvalueDecomposition

/** 
		Determines the eigenvalues and eigenvectors of a real square matrix.
	
 
		If <c>A</c> is symmetric, then <c>A = V * D * V'</c> and <c>A = V * V'</c>
 	where the eigenvalue matrix <c>D</c> is diagonal and the eigenvector matrix <c>V</c> is orthogonal.
 	If <c>A</c> is not symmetric, the eigenvalue matrix <c>D</c> is block diagonal
 	with the real eigenvalues in 1-by-1 blocks and any complex eigenvalues,
 	<c>lambda+i*mu</c>, in 2-by-2 blocks, <c>[lambda, mu; -mu, lambda]</c>.
		The columns of <c>V</c> represent the eigenvectors in the sense that <c>A * V = V * D</c>.
 	The matrix V may be badly conditioned, or even singular, so the validity of the equation
 	<c>A=V*D*inverse(V)</c> depends upon the condition of <c>V</c>.
	
*/
public interface IEigenvalueDecomposition
{
	/** Returns the real parts of the eigenvalues.
	*/
	double[] getRealEigenvalues();

	/** Returns the imaginary parts of the eigenvalues.
	*/
	double[] getImaginaryEigenvalues();

	/** Returns the eigenvector matrix.
	*/
	IMatrix getEigenvectorMatrix();

	/** Returns the block diagonal eigenvalue matrix.
	*/
	IMatrix getDiagonalMatrix();
}