package SignalProcess;


import NumCalc.*;
import ClassicalStat.*;
import java.util.*;

//region Imports


//endregion


/**  Class for convolution calculations 
*/
public class Convolution extends ConvolutionCalc
{

	//region Fields

	private Combinatory comb;
	private Histogram originalHistogram;
	private Histogram histogram;
	private Statistics stat;
	private NormalDist normal;
	private FFT fft;

	private int maxClasses;
	private int maxCombTerms;
	private int maxFftTerms;
	private double mean;
	private double stDev;
	private double n;

	private ArrayList<Double> probs;
	private ArrayList<Double> acumProbs;
	private int lag;

	private int nThreads;

	@FunctionalInterface
	private interface function
	{
		double invoke(int n, double value, boolean acum);
	}

	//endregion

	//region Constructor

	/**  Constructor 
	 @param maxClasses maximum of classes allowed for histogram 
	 @param maxCombTerms maximum of convolution terms for combinatory calculation 
	 @param maxFftTerms maximum of convolution terms for fast fourier calculation 
	 @param nThreads number of threads 
	*/
	public Convolution(int maxClasses, int maxCombTerms, int maxFftTerms, int nThreads)
	{
		this.comb = new Combinatory();
		stat = new Statistics();
		normal = new NormalDist();
		fft = new FFT();

		this.maxClasses = maxClasses;
		this.maxCombTerms = maxCombTerms;
		this.maxFftTerms = maxFftTerms;
		this.nThreads = nThreads;
	}

	//endregion

	//region Properties 

	/**  Maximum number of classes allowed 
	*/
	public final int getMaxClasses()
	{
		return maxClasses;
	}
	public final void setMaxClasses(int value)
	{
		maxClasses = value;
	}

	/**  Maximum number of convolution terms allowed for combinatory 
	*/
	public final int getMaxCombTerms()
	{
		return maxCombTerms;
	}
	public final void setMaxCombTerms(int value)
	{
		maxCombTerms = value;
	}

	/**  Maximum number of convolution terms allowed for fast fourier transformation 
	*/
	public final int getMaxFftTerms()
	{
		return maxFftTerms;
	}
	public final void setMaxFftTerms(int value)
	{
		maxFftTerms = value;
	}

	/**  Histogram for convolution calculation 
	*/
	public final Histogram getHistogram()
	{
		return histogram;
	}

	//endregion

	//region Public Methods

	/**  Load Data interface implementation (needs to set maxClasses first if needed) 
	 @param data
	 @param n
	*/
	public final void LoadData(ArrayList<Double> data, double n)
	{
		Histogram hist = new Histogram(maxClasses);
		hist.LoadData(data);
		this.n = n;
		LoadHistogram(hist, maxClasses);
	}

	/**  Load histogram for calculation 
	 @param hist histogram for calculation 
	 @param n number of convolutions 
	*/
	public final void LoadHistogram(Histogram hist, double n)
	{
		LoadData(hist.GetRawData(), n);
	}

	/**  Load a Histogram for calculation with restriction of classes 
	 @param histogram the histogram 
	 @param maxClasses maximum number of classes allowed (-1 for no restriction) 
	*/
	public final void LoadHistogram(Histogram histogram, int maxClasses)
	{
		this.maxClasses = maxClasses;

		//set properties
		this.originalHistogram = histogram;
		this.mean = originalHistogram.Mean;
		this.stDev = originalHistogram.StDev;
		this.maxClasses = maxClasses;

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


	public final double ProbabilityInt(int n, double value, boolean acum)
	{
		return Probability(n, value, acum);
	}

	/**  ProbabilityComb 
	 @param value value for probability calculation 
	 @param n number of summands 
	 @param acum accumulated probability 
	 @return  probability value 
	*/
	public final double Probability(double n, double value, boolean acum)
	{
		if (histogram.size() == 1)
		{
			double uniqueVal = histogram.GetValues()[0];
			if (acum)
			{
				if (value < uniqueVal)
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
				if (value == uniqueVal)
				{
					return 1;
				}
				else
				{
					return 0;
				}
			}
		}

		if (value >= histogram.Max * n)
		{
			return 1.0;
		}
		if (value <= histogram.Min * n)
		{
			return 0.0;
		}
		int nInt = (int)n;

		//short
		if (n <= maxCombTerms)
		{
			double perm = comb.Permutations(nInt);
			if (nThreads > 1)
			{
				return ProbabilityCombAcumMultithreading(histogram, perm, nInt, value);
			}
			else
			{
				if (acum)
				{
					return ProbabilityCombAcum(histogram, perm, nInt, value);
				}
				else
				{
					return ProbabilityComb(histogram, perm, nInt, value);
				}
			}
		}

		//intermediate
		else if (n <= maxFftTerms)
		{
			if (acum)
			{
				return ProbabilityFftAcum(nInt, value);
			}
			else
			{
				return ProbabilityFft(nInt, value);
			}
		}

		//large
		else
		{
			if (acum)
			{
				return ProbabilityNormalAcum(n, value);
			}
			else
			{
				return 0;
			}
		}
	}

	/**  Quantile 
	 @param n number of summands of the convolution 
	 @param p probability 
	 @return  quantile value 
	*/
	public final double Quantile(double n, double p)
	{
		double quantile = -1;
		int nInt = (int)n;

		if (p == 0.5)
		{
			quantile = histogram.Mean * n;
		}
		else
		{
			//short
			if (nInt <= maxCombTerms)
			{
				if (p < 0.5)
				{
					quantile = MonotoneBisectionRec((int n, double value, boolean acum) -> ProbabilityInt(n, value, acum), histogram.Min * nInt, histogram.Mean * nInt, p, true, 0.01, 1, 100, histogram, nInt);
				}
				else if (p > 0.5)
				{
					quantile = MonotoneBisectionRec((int n, double value, boolean acum) -> ProbabilityInt(n, value, acum), histogram.Mean * nInt, histogram.Max * nInt, p, true, 0.01, 1, 100, histogram, nInt);
				}
			}
			//intermediate
			else if (n <= maxFftTerms)
			{
				quantile = QuantileFft(nInt, p);
			}
			//large
			else
			{
				quantile = QuantileNormal(n, p / 100);
			}
		}
		return Math.round(quantile);
	}

	public final double ProbabilityAcum(double x)
	{
		return Probability(this.n, x, true);
	}

	public final double Quantile(double p)
	{
		return Quantile(this.n, p);
	}


	/**  if calculation is valid 
	 @return  true if it is valid, false if not 
	*/
	public final boolean Valid()
	{
		return true;
	}

	//endregion

	//region Private Methods

	//region Calculation methods

	//region Combinatory (short)

	private double ProbabilityComb(Histogram hist, double perm, int n, double value)
	{
		ArrayList<Integer> keys = hist.GetKeys();
		int key = hist.GetKey(value);
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

	private double ProbabilityCombAcum(Histogram hist, double perm, int n, double value)
	{
		ArrayList<Integer> keys = hist.GetKeys();
		int key = hist.GetKey(value);

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
		for (int i = 0;i < sumsVar.length;i++)
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

	//endregion

	//region Fast Fourier (intermediate)

	private double ProbabilityFft(int n, double value)
	{
		if (probs == null)
		{
			LinearConvolve(histogram, n);
		}
		//int key = (int)(value / histogram.Width- histogram.Min) + lag;
		int key = (int)Math.ceil(((value / histogram.Width - histogram.Min) + histogram.Mean));
		if (key < 0 || key > probs.size() - 1)
		{
			return 0;
		}
		return probs.get(key);
	}

	private double ProbabilityFftAcum(int n, double value)
	{
		if (probs == null)
		{
			LinearConvolve(histogram, n);
		}
		//int key = (int)(value / histogram.Width - histogram.Min) + lag;
		int key = (int)Math.ceil(((value / histogram.Width - histogram.Min + histogram.Mean)));
		if (key < 0 || key > probs.size() - 1)
		{
			return 0;
		}
		return acumProbs.get(key);
	}

	private double QuantileFft(int n, double p)
	{
		if (p > 1)
		{
			throw new RuntimeException(Strings.Error_Probabilities_must_be_between_0_and_1);
		}
		if (probs == null)
		{
			LinearConvolve(histogram, n);
		}

		for (int i = 0;i < acumProbs.size();i++)
		{
			if (acumProbs.get(i).compareTo(p) >= 0)
			{
				if (i > 0 && Math.abs(acumProbs.get(i) - p) > Math.abs(acumProbs.get(i + 1) - p))
				{
					return (int)Math.round(i + 1 - histogram.Mean);
				}
				else
				{
					return (int)Math.round(i - histogram.Mean);
				}
			}
		}
		return -1;
	}

	private void LinearConvolve(Histogram hist, int n)
	{
		ArrayList<Double> completeValues = new ArrayList<Double>();
		for (int i = 0;i < maxClasses;i++)
		{
			if (hist.GetValue(i) > hist.Max)
			{
				break;
			}
			if (hist.Freqs.ContainsKey(i))
			{
				completeValues.add(hist.Freqs[i]);
			}
			else
			{
				completeValues.add(0.0);
			}
		}
		PadRight(completeValues);

		Complex[] comp = new Complex[completeValues.size()];
		for (int i = 0;i < completeValues.size();i++)
		{
			comp[i] = new Complex(completeValues.get(i), 0);
		}
		Complex[] complexRes = fft.LinearConvolve(comp, n);
		ArrayList<Double> realRes = new ArrayList<Double>();
		double tot = 0;
		double val;
		for (int i = 0;i < complexRes.length;i++)
		{
			val = complexRes[i].getReal();
			if (val > 0)
			{
				realRes.add(val);
				tot += val;
			}
			else
			{
				realRes.add(0);
			}

		}
		probs = new ArrayList<Double>();
		acumProbs = new ArrayList<Double>();
		double prob;
		double acum = 0;
		for (int i = 0;i < complexRes.length;i++)
		{
			prob = realRes.get(i) / tot;
			probs.add(realRes.get(i) / tot);
			acum += prob;
			acumProbs.add(acum);
		}
		lag = (int)Math.round(hist.Mean);
	}

	private void CalculateAcum()
	{
		acumProbs = new ArrayList<Double>();
		double acum = 0;
		for (int i = 0;i < probs.size();i++)
		{
			acum += probs.get(i);
			acumProbs.set(i, acum);
		}
	}

	//endregion

	//region Normal (long)

	private double ProbabilityNormalAcum(double n, double value)
	{
		double convMean = n * mean;
		double convStDev = Math.sqrt(n) * stDev;
		if (convStDev == 0)
		{
			if (value < convMean)
			{
				return 0;
			}
			else
			{
				return 1;
			}
		}
		return normal.pNorm(value, convMean, convStDev);
	}

	private double QuantileNormal(double n, double p)
	{
		double convMean = n * mean;
		double convStDev = Math.sqrt(n) * stDev;
		if (convStDev == 0)
		{
			return convMean;
		}
		if (p < 0 || p > 1)
		{
			throw new RuntimeException(Strings.Probability_p_must_be_between_0_and_1);
		}
		return normal.qNorm(p, convMean, convStDev);
	}

	//endregion

	//endregion

	//region Auxiliar Methods

	private double QuantileSec(double min, double max, double acumMin, int n, double p)
	{
		double acum = acumMin;
		double q = 0;
		double acumAnt = -Double.MAX_VALUE;
		for (int i = (int)min + 1;i <= (int)max;i++)
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

	private void PadRight(ArrayList<Double> values)
	{
		int count = values.size();
		int quadPow = 2;
		while (count > quadPow)
		{
			quadPow *= 2;
		}
		if (quadPow == count)
		{
			return;
		}
		int addCount = quadPow - count;
		for (int i = 0;i < addCount;i++)
		{
			values.add(0);
		}
	}

	//endregion

	//region Obsolete Methods

	@Deprecated
	private double ProbabilityCombinations(Histogram hist, int n, double value)
	{
		double perm = comb.Permutations(n);

		ArrayList<Integer> keys = hist.GetKeys();
		ArrayList<int[]> sumsVars = comb.SumsCombinations(keys, n, (int)value, false);
		double convProb = 0;
		double prob;
		for (int[] sumsVar : sumsVars)
		{
			prob = 1;
			for (int i = 0;i < sumsVar.length;i++)
			{
				prob *= hist.ProbabilityByKey(sumsVar[i]);
			}
			convProb += prob;
		}
		return convProb;
	}

	@Deprecated
	private double SimpleConvolutionProbability(Histogram hist, double value)
	{
		int z = (int)value;
		ArrayList<Integer> keys = hist.GetKeys();
		double minValue = keys.get(0);
		double maxValue = keys.get(keys.size() - 1);
		double convProb = 0;
		double Px, Py;
		for (int k = 0;k < keys.size();k++)
		{
			Px = hist.ProbabilityByKey(keys.get(k));
			Py = hist.ProbabilityByKey(z - keys.get(k));
			convProb += Px * Py;
		}
		return convProb;
	}

	@Deprecated
	private double ProbabilityCombValues(Histogram hist, ArrayList<Integer> values, double perm, int n, double value)
	{
		ArrayList<int[]> sumsVars = comb.SumsVariations(values, n, (int)value, false);

		double convProb = 0;
		double prob;
		for (int[] sumsVar : sumsVars)
		{
			prob = ProbabilityOneComb(hist, perm, sumsVar);
			convProb += prob;
		}
		return convProb;
	}

	@Deprecated
	private double ProbabilityCombAcumValues(Histogram hist, ArrayList<Integer> values, double perm, int n, double value)
	{
		HashMap<Integer, ArrayList<int[]>> sumsVarsDict = comb.SumsVariationsUntil(values, n, (int)value, false);

		double convProb = 0;
		double prob;
		for (ArrayList<int[]> sumsVars : sumsVarsDict.values())
		{
			prob = ProbabilityManyCombs(hist, perm, sumsVars);
			convProb += prob;
		}
		return convProb;
	}

	//endregion

	//endregion

	//region Multithreading

	//region Main method

	private double ProbabilityCombAcumMultithreading(Histogram hist, double perm, int n, double value)
	{
		ArrayList<Integer> keys = hist.GetKeys();
		int key = hist.GetKey(value);

		HashMap<Integer, ArrayList<int[]>> sumsVarsDict = comb.SumsVariationsUntil(keys, n, key, false);

		Thread thread;
		MultiCalculator mc;
		ArrayList<Thread> threads = new ArrayList<Thread>();
		ArrayList<MultiCalculator> mcs = new ArrayList<MultiCalculator>();
		for (ArrayList<int[]> sumsVars : sumsVarsDict.values())
		{
			mc = new MultiCalculator(this, sumsVars, hist, perm);
			mcs.add(mc);
			thread = new Thread(new ThreadStart(mc.Calculate), nThreads);
			threads.add(thread);
			thread.start();
		}
		for (Thread th : threads)
		{
			th.join();
		}

		double convProb = 0;
		for (MultiCalculator mC : mcs)
		{
			convProb += mC.getResult();
		}
		return convProb;
	}

	//endregion

	//region Class MultiCalculator

	/**  Inner class for multithreading calculation 
	*/
	public static class MultiCalculator
	{

		//region Fields

		private Convolution conv;
		private double result;
		private ArrayList<int[]> sumsVars;
		private Histogram hist;
		private double perm;

		//endregion

		//region Constructor

		/**  Constructor 
		 @param conv this class, convolution 
		 @param sumsVars sums of variations 
		 @param hist histogram of the distribution 
		 @param perm calculated permutation 
		*/
		public MultiCalculator(Convolution conv, ArrayList<int[]> sumsVars, Histogram hist, double perm)
		{
			this.conv = conv;
			this.sumsVars = sumsVars;
			this.hist = hist;
			this.perm = perm;
		}

		//endregion

		//region Properties

		/**  Result of calculation 
		*/
		public final double getResult()
		{
			return result;
		}

		//endregion

		//region Public Method

		/**  Main calculation method 
		*/
		public final void Calculate()
		{
			result = conv.ProbabilityManyCombs(hist, perm, sumsVars);
		}

		//endregion

		//region Private Methods

		//endregion

	}

	//endregion

	//endregion
}