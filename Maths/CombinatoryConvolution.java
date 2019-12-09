package NumCalc;


//region Imports

import ClassicalStat.*;
import java.util.*;

//endregion



/**  Class for convolution calculation with combinatory 
*/
public class CombinatoryConvolution implements ConvolutionCalc  {

	//region Fields

	private Combinatory comb;
	private int maxClasses = -1;
	private Histogram originalHistogram;
	private Histogram histogram;
	private double mean;
	private double stDev;
	private double n;
	private boolean acum;
	private double maxError;
	@FunctionalInterface
	private interface function
	{
		double invoke(int n, double value, boolean acum);
	}

	//endregion

	//region Constructor

	/**  Constructor 
	 @param maxClasses maximum of classes of the histogram (scaling) 
	*/
	public CombinatoryConvolution(int maxClasses, double maxError)
	{
		this.maxClasses = maxClasses;
		this.comb = new Combinatory();
		this.acum = true;
		this.maxError = maxError;
	}

	//endregion

	//region Public Methods

	/**  Load data for calculation 
	 @param data data from time series 
	 @param n number (real) of convolutions 
	*/
	public final void LoadData(ArrayList<Double> data, double n)
	{
		Histogram hist = new Histogram(maxClasses);
		hist.LoadData(data);
		this.n = n;
		LoadHistogram(hist, n);
	}

	/**  Load a Histogram for calculation with restriction of classes 
	 @param histogram the histogram 
	 @param n number (real) of convolutions 
	*/
	public final void LoadHistogram(Histogram histogram, double n)
	{

		//set properties
		this.originalHistogram = histogram;
		this.mean = originalHistogram.Mean;
		this.stDev = originalHistogram.StDev;
		this.n = n;

		//scale classes
		if (maxClasses > 0 && maxClasses < histogram.Freqs.size())
		{
			this.histogram = new Histogram(maxClasses);
			this.histogram.LoadHist(originalHistogram);
		}
		else
		{
			this.histogram = this.originalHistogram;
		}
	}

	/**  Acumulated probability for a certain value 
	 @param x the value 
	 @return  the probability 
	*/
	public final double ProbabilityAcum(double x)
	{
		if (histogram.size() == 1)
		{
			double uniqueVal = histogram.GetValues()[0];
			if (acum)
			{
				if (x < uniqueVal)
				{
					return 0;
				}
				else
				{
					return 1;
				}
			}
			else
			{
				if (x == uniqueVal)
				{
					return 1;
				}
				else
				{
					return 0;
				}
			}
		}

		if (x >= histogram.Max * n)
		{
			return 1.0;
		}
		if (x <= histogram.Min * n)
		{
			return 0.0;
		}
		int nInt = (int)n;

		double perm = comb.Permutations(nInt);
		if (acum)
		{
			return ProbabilityAcum(histogram, perm, nInt, x);
		}
		else
		{
			return Probability(histogram, perm, nInt, x);
		}


	}

	/**  Quantile for a certain probability 
	 @param p the probability 
	 @return  the value 
	*/
	public final double Quantile(double p)
	{
   double quantile = -1;
		int nInt = (int)n;
		double min = histogram.Min * nInt;
		if (min < 0)
		{
			min = 0;
		}
		quantile = MonotoneBisectionRec((int n, double value, boolean acum) -> Probability(n, value, acum), min, histogram.Max * nInt, p, true, 0.01, 1, 100, histogram, nInt);
		return Math.round(quantile);
	}

	public final double Probability(int n, double x, boolean acum)
	{
		return ProbabilityAcum(x);
	}

	/**  if calculation is valid 
	 @return  true if it is valid, false if not 
	*/
	public final boolean Valid()
	{
		double quantMean = Quantile(0.5);
		return (Math.abs(quantMean - mean * n) / (mean * n) <= maxError);
	}

	//endregion

	//region Private Methods

	private double Probability(Histogram hist, double perm, int n, double value)
	{
		ArrayList<Integer> keys = hist.GetKeys();
		int key = hist.GetConvKey(value, n);
		ArrayList<int[]> sumsVars = comb.SumsVariations(keys, n, key, false);

		double convProb = 0;
		double prob;
		for (int[] sumsVar : sumsVars)
		{

			prob = ProbabilityOneComb(hist, perm, sumsVar);
			convProb += prob;
		}
		return convProb;
	}

	private double ProbabilityAcum(Histogram hist, double perm, int n, double value)
	{
		ArrayList<Integer> keys = hist.GetKeys();
		int key = hist.GetConvKey(value, n);

		HashMap<Integer, ArrayList<int[]>> sumsVarsDict = comb.SumsVariationsUntil(keys, n, key, false);

		double convProb = 0;
		double prob;
		for (ArrayList<int[]> sumsVars : sumsVarsDict.values())
		{
			prob = ProbabilityManyCombs(hist, perm, sumsVars);
			convProb += prob;
		}
		return convProb;
	}

	private double ProbabilityOneComb(Histogram hist, double perm, int[] sumsVar)
	{
		double prob = 1;
		int rep = 1;
		double permRep = 1;
		for (int i = 0; i < sumsVar.length; i++)
		{
			if (i > 0 && sumsVar[i] == sumsVar[i - 1])
			{
				rep++;
			}
			if (rep > 1 && (i == sumsVar.length - 1 || sumsVar[i] != sumsVar[i - 1]))
			{
				permRep *= comb.Permutations(rep);
				rep = 1;
			}
			prob *= hist.ProbabilityByKey(sumsVar[i]);
		}
		return prob * (perm / permRep);
	}

	private double ProbabilityManyCombs(Histogram hist, double perm, ArrayList<int[]> sumsVars)
	{
		double oneProb = 0;
		double manyProbs = 0;
		for (int[] sumsVar : sumsVars)
		{
			oneProb = ProbabilityOneComb(hist, perm, sumsVar);
			manyProbs += oneProb;
		}
		return manyProbs;
	}

	private double QuantileSec(double min, double max, double acumMin, int n, double p)
	{
		double acum = acumMin;
		double q = 0;
		double acumAnt = -Double.MAX_VALUE;
		for (int i = (int)min + 1; i <= (int)max; i++)
		{
			q = i;
			acumAnt = acum;
			acum += Probability(n, q, false);
			if (acum > p)
			{
				break;
			}
		}
		if (Math.abs(acumAnt - p) < Math.abs(acum - p))
		{
			return q - 1;
		}
		else
		{
			return q;
		}
	}

	//TODO: Refactoring (this is almost a clone of MonotoneBisectionRec of NumCalc (cannot use the other due to circular reference)
	private double MonotoneBisectionRec(function f, double min, double max, double val, boolean ascendant, double eps, int it, int maxIt, Histogram hist, int n)
	{
		it++;
		double inter = min + (max - min) / 2.0;
		if (it > maxIt || inter == min || inter == max)
		{
			return inter;
		}
		double fInter = f.invoke(n, inter, true);
		double diff = val - fInter;
		if (Math.abs(diff) < eps)
		{
			return inter;
		}
		if (diff > 0 && diff < 0.05)
		{
			return QuantileSec(inter, max, fInter, n, val);
		}
		if (ascendant)
		{
			if (val < fInter)
			{
				return MonotoneBisectionRec(f, min, inter, val, ascendant, eps, it, maxIt, hist, n);
			}
			if (val > fInter)
			{
				return MonotoneBisectionRec(f, inter, max, val, ascendant, eps, it, maxIt, hist, n);
			}
		}
		else
		{
			if (val > fInter)
			{
				return MonotoneBisectionRec(f, min, inter, val, ascendant, eps, it, maxIt, hist, n);
			}
			if (val < fInter)
			{
				return MonotoneBisectionRec(f, inter, max, val, ascendant, eps, it, maxIt, hist, n);
			}
		}
		return inter;
	}


	//endregion

}