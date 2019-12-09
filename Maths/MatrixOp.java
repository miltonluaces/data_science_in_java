package NumCalc;

import NumCalc.Properties.*;

/**  Clase de operaciones para matrices 
*/
public class MatrixOp
{

	//region Funciones Publicas

	/**  Numero de filas de una matriz 
	 @param mat the matrix 
	 @return  number of rows 
	*/
	public static int GetRows(double[][] mat)
	{
		return mat.length;
	}

	/**  Numero de columnas de una matriz 
	 @param mat the matrix 
	 @return  number of columns 
	*/
	public static int GetCols(double[][] mat)
	{
		return mat[0].length;
	}

	/**  Suma de matrices 
	 @param A first matrix 
	 @param B second matrix 
	 @return  result matrix 
	*/
	public static double[][] Addition(double[][] A, double[][] B)
	{
		double[][] sum = new double[A.length][A[0].length];
		for (int i = 0;i < A.length;i++)
		{
			for (int j = 0;j < A[0].length;j++)
			{
				sum[i][j] = A[i][j] + B[i][j];
			}
		}
		return sum;
	}

	/**  Resta de matrices 
	 @param A first matrix 
	 @param B second matrix 
	 @return  result matrix 
	*/
	public static double[][] Substraction(double[][] A, double[][] B)
	{
		double[][] subs = new double[A.length][A[0].length];
		for (int i = 0;i < A.length;i++)
		{
			for (int j = 0;j < A[0].length;j++)
			{
				subs[i][j] = A[i][j] - B[i][j];
			}
		}
		return subs;
	}

	/**  Multiplicacion de matrices 
	 @param A first matrix 
	 @param B second matrix 
	 @return  result matrix 
	*/
	public static double[][] Multiply(double[][] A, double[][] B)
	{
		//if(IsStrassenLike(A) && IsStrassenLike(B)) { return MultiplyStrassen(A, B); }
		int fil = A.length;
		int col = B[0].length;
		int p = A[0].length;

		double[][] prod = new double[fil][col];
		for (int i = 0;i < fil;i++)
		{
			for (int j = 0;j < col;j++)
			{
				prod[i][j] = 0.0;
				for (int k = 0;k < p;k++)
				{
					prod[i][j] += A[i][k] * B[k][j];
				}
			}
		}
		return prod;
	}

	/**  Multiplicacion de triangulares 
	 @param A first matrix 
	 @param B second matrix 
	 @return  result matrix 
	*/
	public static double[][] MultiplyT(double[][] A, double[][] B)
	{
		//int row = A.GetLength(0);
		//int col = B.GetLength(1);
		//double[,] prod = new double[row, col];
		//prod = MatrixOp.Multiply(A, B);
		//MatrixOp.MakeSimetric(prod);
		//return (prod);

		int p = A[0].length;
		int fil = A.length;
		int col = B[0].length;
		double[][] prod = new double[fil][col];
		for (int i = 0;i < fil;i++)
		{
			for (int j = i;j < col;j++)
			{
				prod[i][j] = 0.0;
				for (int k = 0;k < p;k++)
				{
					prod[i][j] += A[j][k] * B[k][i];
				}
				prod[j][i] = prod[i][j];
			}
		}
		return prod;
	}

	/**  Matriz Traspuesta 
	 @param mat the matrix 
	 @return  the trasposed matrix 
	*/
	public static double[][] Transpose(double[][] mat)
	{
		double[][] trasp = new double[mat[0].length][mat.length];
		for (int i = 0;i < mat.length;i++)
		{
			for (int j = 0;j < mat[0].length;j++)
			{
				trasp[j][i] = mat[i][j];
			}
		}
		return trasp;
	}

	/**  Multiplicacion por su traspuesta 
	 @param A the matrix 
	 @return  the multiplied trasposed matrix 
	*/
	public static double[][] MultiplyTranspose(double[][] A)
	{
		double[][] tr = Transpose(A);
		double[][] prod = Multiply(A, tr);
		return prod;
	}

	/**  Multiplicacion por un escalar 
	 @param mat the matrix 
	 @param val the scalar 
	 @return  the product matrix 
	*/
	public static double[][] Multiply(double[][] mat, double val)
	{
		double[][] prod = new double[mat.length][mat[0].length];
		for (int i = 0;i < mat.length;i++)
		{
			for (int j = 0;j < mat[0].length;j++)
			{
				prod[i][j] = mat[i][j] * val;
			}
		}
		return prod;
	}

	/**  Hacer simetrica una matriz 
	 @param mat the matrix 
	*/
	public static void MakeSimetric(double[][] mat)
	{
		int row = mat.length;
		int col = mat[0].length;
		for (int c = 0;c < col;c++)
		{
			for (int r = c;r < row;r++)
			{
				mat[c][r] = mat[r][c];
			}
		}
	}

	/**  Clone 
	 @param mat the matrix 
	 @return  the cloned matrix 
	*/
	public static double[][] Clone(double[][] mat)
	{
		double[][] clone = new double[mat.length][mat[0].length];
		for (int i = 0;i < mat.length;i++)
		{
			for (int j = 0;j < mat[0].length;j++)
			{
				clone[i][j] = mat[i][j];
			}
		}
		return clone;
	}

	/**  Paso a escalar de matrices de 1 x 1 
	 @param mat the matrix 
	 @return  the scalar 
	*/
	public static double ToScalar(double[][] mat)
	{
		if (mat.length != 1 || mat[0].length != 1)
		{
			throw new RuntimeException(Strings.No_es_de_1x1);
		}
		return mat[0][0];
	}

	/**  Si todos los elementos son cero 
	 @param mat the matrix 
	 @return  if all are zeros 
	*/
	public static boolean IsZero(double[][] mat)
	{
		double sum = 0.0;
		for (int i = 0;i < mat.length;i++)
		{
			for (int j = 0;j < mat[0].length;j++)
			{
				sum += mat[i][j];
			}
		}
		return sum == 0;
	}

	//endregion

	//region Funciones Privadas

	/**  Multiplicacion performante de cuadradas. Algoritmo de Strassen 
	 @param A first matrix 
	 @param B second matrix 
	 @return  result matrix 
	*/
	public static double[][] MultiplyStrassen(double[][] A, double[][] B)
	{
		double[][] prod = new double[A.length][B[0].length];

		//caso base
		if (A.length == 2)
		{
			prod = Multiply(A, B);
		}
			//caso recursivo
		else
		{
			double[][] A11 = GetCuadrante(A, 1);
			double[][] B11 = GetCuadrante(B, 1);
			double[][] A12 = GetCuadrante(A, 2);
			double[][] B12 = GetCuadrante(B, 2);
			double[][] A21 = GetCuadrante(A, 3);
			double[][] B21 = GetCuadrante(B, 3);
			double[][] A22 = GetCuadrante(A, 4);
			double[][] B22 = GetCuadrante(B, 4);

			double[][] P = MultiplyStrassen(Addition(A11, A22), Addition(B11, B22));
			double[][] Q = MultiplyStrassen(Addition(A21, A22), B11);
			double[][] R = MultiplyStrassen(A11, Substraction(B12, B22));
			double[][] S = MultiplyStrassen(A22, Substraction(B21, B11));
			double[][] T = MultiplyStrassen(Addition(A11, A12), B22);
			double[][] U = MultiplyStrassen(Substraction(A21, A11), Addition(B11, B12));
			double[][] V = MultiplyStrassen(Substraction(A12, A22), Addition(B21, B22));

			double[][] C11 = Substraction(Addition(Addition(P, S), V), T);
			double[][] C12 = Addition(R, T);
			double[][] C21 = Addition(Q, S);
			double[][] C22 = Substraction(Addition(Addition(P, R), U), Q);

			prod = Concat(C11, C12, C21, C22);
		}
		return prod;
	}

	private static double[][] GetCuadrante(double[][] mat, int pos)
	{
		int lado2 = mat.length;
		int lado = lado2 / 2;
		double[][] cuad = new double[lado][lado];

		switch (pos)
		{
			case 1:
				for (int i = 0;i < lado;i++)
				{
					for (int j = 0;j < lado;j++)
					{
						cuad[i][j] = mat[i][j];
					}
				}
				break;
			case 2:
				for (int i = 0;i < lado;i++)
				{
					for (int j = lado;j < lado2;j++)
					{
						cuad[i][j - lado] = mat[i][j];
					}
				}
				break;
			case 3:
				for (int i = lado;i < lado2;i++)
				{
					for (int j = 0;j < lado;j++)
					{
						cuad[i - lado][j] = mat[i][j];
					}
				}
				break;
			case 4:
				for (int i = lado;i < lado2;i++)
				{
					for (int j = lado;j < lado2;j++)
					{
						cuad[i - lado][j - lado] = mat[i][j];
					}
				}
				break;
		}
		return cuad;
	}

	private static double[][] Concat(double[][] A, double[][] B, double[][] C, double[][] D)
	{
		int lado = A.length;
		double[][] union = new double[lado * 2][lado * 2];

		for (int i = 0;i < lado;i++)
		{
			for (int j = 0;j < lado;j++)
			{
				union[i][j] = A[i][j];
				union[i][j + lado] = B[i][j];
				union[i + lado][j] = C[i][j];
				union[i + lado][j + lado] = D[i][j];
			}
		}
		return union;
	}

	private static boolean IsStrassenLike(double[][] mat)
	{
		int i = 1;
		while (i <= 30)
		{
			i *= 2;
			if (mat.length == mat[0].length && mat.length == i)
			{
				return true;
			}
		}
		return false;
	}


	//endregion

}