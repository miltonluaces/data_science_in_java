package NumCalc;

import NumCalc.Properties.*;
import ClassicalStat.*;

//region Imports


//endregion


//region Interfaces 

//region Interface IMatrix

/** Matrix provides the fundamental operations of numerical linear algebra.
*/
public interface IMatrix
{
	/** Returns the number of columns.
	*/
	int getRows();

	/** Returns the number of columns.
	*/
	int getColumns();

	/** Access the value at the given location. First index is Row, second index is column.
	*/
	double getItem(int i, int j);
	void setItem(int i, int j, double value);

	/** Sets all cells of a matrix to value v.
	 @param v Value to assign to all cells of the matrix.
	 @return  if result is ok 
	*/
	boolean Set_Matrix(double v);

	/** Returns a sub matrix extracted from the current matrix.
	 @param startRow Start row index.
	 @param endRow End row index;
	 @param startColumn Start column index;
	 @param endColumn End column index;
	 @return  the calculated matrix 
	*/
	IMatrix Submatrix(int startRow, int endRow, int startColumn, int endColumn);

	/** Returns a sub matrix extracted from the current matrix.
	 @param r Array of row indices;
	 @param c Array of row indices;
	 @return  the calculated matrix 
	*/
	IMatrix Submatrix(int[] r, int[] c);

	/** Returns a sub matrix extracted from the current matrix.
	 @param startRow Starttial row index.
	 @param endRow End row index.
	 @param c Array of row indices.
	 @return  the calculated matrix 
	*/
	IMatrix Submatrix(int startRow, int endRow, int[] c);

	/** Returns a sub matrix extracted from the current matrix.
	 @param r Array of row indices.
	 @param startColumn Start column index.
	 @param endColumn End column index.
	 @return  the calculated matrix 
	*/
	IMatrix Submatrix(int[] r, int startColumn, int endColumn);

	/** Creates a copy of the matrix.
	 @return  the calculated matrix 
	*/
	IMatrix Clone();

	/** Returns the transposed matrix.
	 @return  the calculated matrix 
	*/
	IMatrix Transpose();

	/** Inverse of the matrix if matrix is square, pseudoinverse otherwise.
	*/
	IMatrix getInverse();

	/** Determinant if matrix is square.
	*/
	double getDeterminant();

	/** Returns the One Norm for the matrix.
	 <value>The maximum column sum.</value>
	*/
	double getNorm1();

	/** Returns the Infinity Norm for the matrix.
	 <value>The maximum row sum.</value>
	*/
	double getInfinityNorm();

	/** Returns the Frobenius Norm for the matrix.
	 <value>The square root of sum of squares of all elements.</value>
	*/
	double getFrobeniusNorm();

	/** Return <see langword="true"/> if the matrix is a square matrix.
	*/
	boolean getIsSquare();

	/** Returns <see langword="true"/> if the matrix is symmetric.
	*/
	boolean getIsSymmetric();

	/** Return <see langword="true"/> if the matrix has all its value to 0.0
	*/
	boolean getIsZero();

	/** Returns the trace of the matrix.
	 @return Sum of the diagonal elements.
	*/
	double getTrace();

	/** Matrix addition.
	 @param B the matrix to add 
	 @return  result matrix 
	*/
	IMatrix Addition(IMatrix B);

	/** Method added by ECC. Convert from Matrix to scalar.
	 @return  the scalar 
	*/
	double Toscalar();

	/** Method added by ECC. Convert from Matrix to scalar.
	 @return  result matrix 
	 @param precision number of decimals 
	*/
	IMatrix Round(int precision);

	/** Matrix-matrix multiplication.
	 @param B the matrix to multiply 
	 @return  result matrix 
	*/
	IMatrix Multiply(IMatrix B);

	/** 
	 Matrix-matrix triangular multiplication.
	 Just Multiplies to get the left down diagonal and then copies to the
	 right upper diagonal.
	 
	 @param B the matrix to multiply 
	 @return  result matrix 
	*/
	IMatrix MultiplyT(IMatrix B);

	/** 
	 Copies the left botton corner (triangle) of the matrix to the right upper corner.
	*/
	void CopyToRightUpperCorner();

	/** 
	 Makes this matrix simetric by copying the right upper diagonal into.
	 the left down diagonal
	 
	 @return  result matrix 
	*/
	IMatrix MakeSimetric();

	/** Multiply Matrix by it transpose.
	 @param reverse If true, Multiply Transpose by Matrix. If False, Matrix by its Transpose
	 @param precision Precision to round numbers
	 @return  result matrix 
	*/
	IMatrix MultiplyTranspose(boolean reverse, int precision);

	/** Matrix-scalar multiplication.
	 @param s the scalar to multiply 
	 @return  result matrix 
	*/
	IMatrix Multiply(double s);

	/** Matrix subtraction.
	 @param B the matrix for the operation 
	 @return  result matrix 
	*/
	IMatrix Subtraction(IMatrix B);

	/** Returns the LHS solution vetor if the matrix is square or the least squares solution otherwise.
	 @param rhs Right hand side matrix with as many rows as <c>A</c> and any number of columns.
	 @return  result matrix 
	*/
	IMatrix Solve(IMatrix rhs);

	/** Returns matrix of doubles.
	 @return  result matrix 
	*/
	double[][] ConvertToDouble();

	/** Returns the cholesky decomposition for this matrix.
	 @return  result descomposition 
	*/
	ICholeskyDecomposition GetCholeskyDecomposition();

	/** Returns the LU decomposition for this matrix.
	 @return  result descomposition 
	*/
	ILuDecomposition GetLuDecomposition();

	/** Returns the singular value decomposition for this matrix.
	 @return  result descomposition 
	*/
	ISingularValueDecomposition GetSingularValueDecomposition();

	/** Returns the QR decomposition for this matrix.
	 @return  result descomposition 
	*/
	IQrDecomposition GetQrDecomposition();

	/** Returns the eigenvalue decomposition for this matrix.
	 @return  result descomposition 
	*/
	IEigenvalueDecomposition GetEigenvalueDecomposition();
}