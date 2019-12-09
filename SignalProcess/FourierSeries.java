package SignalProcess;

import java.util.*;

//region Imports


//endregion


/**  Class FourierSeries for fourier discrete trigonometric polynomial regression 
*/
public class FourierSeries extends SignalProc
{

	//region Fields

	private ArrayList<Double> X;
	private ArrayList<Double> Y;
	private ArrayList<Double> aCoeffs;
	private ArrayList<Double> bCoeffs;
	private int n;
	private int m;

	//endregion

	//region Constructor

	/**  Constructor 
	*/
	public FourierSeries()
	{
		aCoeffs = new ArrayList<Double>();
		bCoeffs = new ArrayList<Double>();
	}

	//region Properties

	/**  Number of coefficients 
	*/
	public final int getN()
	{
		return n;
	}

	/**  A fourier coefficents 
	*/
	public final ArrayList<Double> getACoeffs()
	{
		return aCoeffs;
	}

	/**  B fourier coefficents 
	*/
	public final ArrayList<Double> getBCoeffs()
	{
		return bCoeffs;
	}


	//endregion

	//endregion

	//region Public Methods

	//region Load Data

	/**  Load Data for fourier trigonometric polynomial regression 
	 @param X uniform distributed independent values 
	 @param Y dependent values 
	*/
	public final void LoadData(ArrayList<Double> X, ArrayList<Double> Y)
	{
		if (X.size() != Y.size())
		{
			throw new RuntimeException(Strings.X_must_be_of_same_size_than_Y);
		}
		if (X.size() % 2 != 0)
		{
			throw new RuntimeException(Strings.X_Y_count_must_be_even);
		}
		this.X = X;
		this.Y = Y;
		this.m = X.size() / 2;
	}

	/**  Load data for calculation 
	 @param Y list of data 
	*/
	public final void LoadData(ArrayList<Double> Y)
	{
		if (Y.size() % 2 != 0)
		{
			Y.add(0, Y.get(0));
		}
		this.m = Y.size() / 2;
		this.X = new ArrayList<Double>();
		for (int j = 0;j < 2 * m;j++)
		{
			X.add(-Math.PI + ((double)j / (double)m) * Math.PI);
		}
		this.Y = Y;
	}

	//endregion

	//region Coefficents

	/**  Calculate n fourier coefficients ak, bk 
			  2m-1
	 ak = 1/m ∑ yj cos kxj
			  j=0
	
			  2m-1
	 bk = 1/m ∑ yj sin kxj
			  j=0
	 @param n number of coefficients 
	*/
	public final void CalculateCoeffs(int n)
	{
		this.n = n;
		this.aCoeffs.clear();
		this.bCoeffs.clear();
		double a, b;

		for (int k = 0;k <= n;k++)
		{
			a = 0.0;
			for (int j = 0;j < 2 * m;j++)
			{
				a += Y.get(j) * Math.cos(k * X.get(j));
			}
			a = a / (double)m;
			aCoeffs.add(a);
		}

		bCoeffs.add(0.00); //no bCoeffs in zero pos
		for (int k = 1;k < n;k++)
		{
			b = 0.0;
			for (int j = 0;j < 2 * m;j++)
			{
				b += Y.get(j) * Math.sin(k * X.get(j));
			}
			b = b / (double)m;
			bCoeffs.add(b);
		}
	}

	/**  Set directly coefficents 
	 @param aCoeffs cos coefficents 
	 @param bCoeffs sin coefficents 
	*/
	public final void SetCoeffs(ArrayList<Double> aCoeffs, ArrayList<Double> bCoeffs)
	{
		aCoeffs.clear();
		bCoeffs.clear();
		for (int k = 0;k < aCoeffs.size();k++)
		{
			aCoeffs.add((Math.pow(-1, k) / (double)m) * aCoeffs.get(k));
			bCoeffs.add((Math.pow(-1, k) / (double)m) * bCoeffs.get(k));
		}
		this.aCoeffs = aCoeffs;
		this.bCoeffs = bCoeffs;
		this.n = aCoeffs.size() - 1;
	}

	//endregion

	//region Calculate Result

	/**  Fourier Regression Function 
							  n-1  
	s(x) = a0/2 + an cos nx + ∑ (ak cos kx + bk sen kx)
							  df=1      
	 @param x independent variable 
	 @return  y dependent variable 
	*/
	public final double Function(double x)
	{
		double res = aCoeffs.get(0) / 2.0 + aCoeffs.get(n) * Math.cos(n * x);
		for (int k = 1;k < n;k++)
		{
			res += (aCoeffs.get(k) * Math.cos(k * x) + bCoeffs.get(k) * Math.sin(k * x));
		}
		return res;
	}

	/**  Calculate whole series values 
	 @return  calculated serie values 
	*/
	public final ArrayList<Double> CalculateSerie()
	{
		ArrayList<Double> serie = new ArrayList<Double>();
		for (double x : X)
		{
			serie.add(Function(x));
		}
		return serie;
	}

	//endregion

	//endregion

}