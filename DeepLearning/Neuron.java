package DeepLearning;

import java.util.*;

//region Imports


//endregion


public class Neuron implements Serializable
{

	//region Fields

	private static Random rand = new Random();

	public int Rows;
	public int Cols;
	public double[] W;
	public double[] Dw;
	public double[] Cache;

	//endregion

	//region Constructors

	public Neuron(int rows, int cols)
	{
		this.Rows = rows;
		this.Cols = cols;
		this.W = new double[Rows * Cols];
		this.Dw = new double[Rows * Cols];
		this.Cache = new double[Rows * Cols];
	}

	public Neuron(int dim)
	{
		this(dim, 1);
	}

	public Neuron(double[] vector)
	{
		this(vector.length, 1);
		this.W = vector;
	}

	//endregion

	//region Public Methods

	//region General

	public final Neuron Clone()
	{
		Neuron res = new Neuron(Rows, Cols);
		for (int i = 0; i < W.length; i++)
		{
			res.W[i] = W[i];
			res.Dw[i] = Dw[i];
			res.Cache[i] = Cache[i];
		}
		return res;
	}

	@Override
	public String toString()
	{
		String res = "";
		for (int r = 0; r < Rows; r++)
		{
			for (int c = 0; c < Cols; c++)
			{
//C# TO JAVA CONVERTER TODO TASK: The '0:N5' format specifier is not converted to Java:
				res += String.format("{0:N5}", GetW(r, c)) + "\t";
			};
			res += "\n";
		}
		return res;
	}

	//endregion

	//region SetGetters

	private int GetIndex(int row, int col)
	{
		int ix = row * this.Cols + col;
		return ix;
	}

	private double GetW(int row, int col)
	{
		return W[GetIndex(row, col)];
	}

	private void SetW(int row, int col, double val)
	{
		W[GetIndex(row, col)] = val;
	}

	//endregion

	//region Reset

	public final void ResetDw()
	{
		for (int i = 0; i < Dw.length; i++)
		{
			Dw[i] = 0;
		}
	}

	public final void ResetCache()
	{
		for (int i = 0; i < Cache.length; i++)
		{
			Cache[i] = 0;
		}
	}

	//endregion

	//region Operations

	public static Neuron Transpose(Neuron m)
	{
		Neuron res = new Neuron(m.Cols, m.Rows);
		for (int r = 0; r < m.Rows; r++)
		{
			for (int c = 0; c < m.Cols; c++)
			{
				res.SetW(c, r, m.GetW(r, c));
			}
		}
		return res;
	}

	public static Neuron Random(int rows, int cols, double stDev)
	{
		Neuron res = new Neuron(rows, cols);
		for (int i = 0; i < res.W.length; i++)
		{
			res.W[i] = rand.nextDouble() * stDev;
		}
		return res;
	}

	public static Neuron Ident(int dim)
	{
		Neuron res = new Neuron(dim, dim);
		for (int i = 0; i < dim; i++)
		{
			res.SetW(i, i, 1.0);
		}
		return res;
	}

	public static Neuron Rep(int rows, int cols, double val)
	{
		Neuron res = new Neuron(rows, cols);
		for (int i = 0; i < res.W.length; i++)
		{
			res.W[i] = val;
		}
		return res;
	}

	//endregion


	//endregion
}