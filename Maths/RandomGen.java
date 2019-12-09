package NumCalc;

//region Imports

import ClassicalStat.*;
import java.util.*;

//endregion


/** 
 RandomProvider.  Provides random numbers of all data types in specified ranges.  It also contains a couple of methods
 from Normally (Gaussian) distributed random numbers and Exponentially distributed random numbers.
*/
public class RandomGen
{

	//region Fields

	private Random rand;
	private double storedUniformDeviate;
	private boolean storedUniformDeviateIsGood = false;
	private Functions func;

	private static int seed = 0;

	//endregion

	//region Constructor

	/**  Constructor 
	*/
	public RandomGen()
	{
		func = new Functions();
		Reset();
	}

	/**  Reset seed 
	*/
	public final void Reset() {
		rand = new Random(seed == 0 ? System.currentTimeMillis() : seed);
	}

	/**  Force seed. 
	 @param seed seed for next random numbers set.
	*/
	public static void SetSeed(int newSeed)
	{
		seed = newSeed;
	}

	//endregion

	//region Distributions

	//region Uniform

	/**  Devuelve un double rand�mico en el rango [0,1) (cero inc) 
	 @return  result 
	*/
	public final double NextDouble()
	{
		double nextDouble = 0.0;
		synchronized (rand)
		{
			nextDouble = rand.nextDouble();
		}
		return nextDouble;
	}

	/**  Devuelve un booleano rand�mico 
	 @return  result 
	*/
	public final boolean NextBoolean()
	{
		double nextBoolean = 0.0;
		synchronized (rand)
		{
			nextBoolean = rand.nextInt(2);
		}
		return nextBoolean != 0;
	}

	/**  Devuelve un int en el rango [min, max) 
	 @param min minimum value 
	 @param max maximum value 
	 @return  result 
	*/
	public final int NextInt(int min, int max)
	{
		if (max == min)
		{
			return min;
		}
		else if (max < min)
		{
			throw new IllegalArgumentException("Max_must_be_greater_than_min");
		}
		double nextDouble = 0.0;
		synchronized (rand)
		{
			nextDouble = rand.nextDouble();
		}
		int randInt = (int)java.lang.Math.round(min + nextDouble * (max - min));
		return randInt;
	}

	/**  Devuelve un double en el rango [min, max) 
	 @param min minimum value 
	 @param max maximum value 
	 @return  result 
	*/
	public final double NextDouble(double min, double max)
	{
		if (max <= min)
		{
			throw new IllegalArgumentException("Max_must_be_greater_than_min");
		}
		double nextDouble = 0.0;
		synchronized (rand)
		{
			nextDouble = rand.nextDouble();
		}
		double randDbl = min + nextDouble * (max - min);
		return randDbl;
	}

	/**  Devuelve un double en el rango [0,1) con distribucion Uniforme 
	 @return  result 
	*/
	public final double NextUniform()
	{
		return NextDouble();
	}

	/**  Devuelve un double en el rango dado con distribucion Uniforme 
	 @param min start of interval 
	 @param max end of interval 
	 @return  result 
	*/
	public final double NextUniform(double min, double max)
	{
		return NextDouble(min, max);
	}

	/**  Generar array de randomicos 
	 @param cantElementos number of elements 
	 @return  random array 
	*/
	public final double[] GenerarArrayRandom(int cantElementos)
	{
		double[] array = new double[cantElementos];
		double nextDouble = 0.0;
		synchronized (rand)
		{
			nextDouble = rand.nextInt();
		}
		for (int i = 0;i < cantElementos;i++)
		{
			array[i] = nextDouble;
		}
		return array;
	}

	//endregion

	//region Normal

	/**  Devuelve variables tipificadas N(0,1) (sesgos) 
	 @return  next normal random value 
	*/
	public final double NextNormal()
	{
		// basado en algoritmo de Numerical Recipes
		if (storedUniformDeviateIsGood)
		{
			storedUniformDeviateIsGood = false;
			return storedUniformDeviate;
		}
		else
		{
			double rsq = 0.0;
			double v1 = 0.0, v2 = 0.0, fac = 0.0;
			while (rsq == 0.0 || rsq >= 1.0)
			{
				v1 = NextDouble() * 2.0 - 1.0;
				v2 = NextDouble() * 2.0 - 1.0;
				rsq = Math.pow(v1, 2) + Math.pow(v2, 2);
			}
			fac = Math.sqrt(Math.log(rsq) / rsq * -2.0);
			storedUniformDeviate = v1 * fac;
			storedUniformDeviateIsGood = true;
			return v2 * fac;
		}
	}

	/**  Devuelve variables N(m,s) (sesgos) 
	 @param m mean 
	 @param s standard deviation 
	 @return  next normal random value 
	*/
	public final double NextNormal(double m, double s)
	{
		double z = NextNormal();
		double x = s * z + m;
		return x;
	}


	//endregion

	//region Exponential

	/**  Devuelve sesgos rand�micos positivos con media = 1, con distribucion Exponencial 
	 @return  next exponential random value 
	*/
	public final double NextExponential()
	{
		double dum = 0.0;
		while (dum == 0.0)
		{
			dum = NextUniform();
		}
		return -Math.log(dum);
	}

	/**  Devuelve sesgos rand�micos positivos con media = m, con distribucion Exponencial 
	 @param m media value 
	 @return  next exponential random value 
	*/
	public final double NextExponential(double m)
	{
		return NextExponential() + m;
	}

	//endregion

	//region Empirical

	/**  Get a random value with the histogram distribution 
	 @param hist the histogram that contains the distribution 
	 @return  the random value 
	*/
	public final double GetRandomValue(Histogram hist)
	{

		double val = NextDouble(0, hist.TotFreqs);
		int valIndex = BinarySearch(hist.AccumFreqValues, val);
		int key = hist.AccumFreqKeys[valIndex];
		return hist.GetValue(key);
	}

	/**  Get a secquence of n random values with the histogram distribution 
	 @param hist the histogram that contains the distribution 
	 @param n number of values requiered 
	 @return  the random values 
	*/
	public final ArrayList<Double> GetRandomValues(Histogram hist, int n)
	{
		ArrayList<Double> vals = new ArrayList<Double>();
		for (int i = 0;i < n;i++)
		{
			vals.add(GetRandomValue(hist));
		}
		return vals;
	}

	/**  Get a secquence of n random values with the histogram distribution with a certain likelyhood 
	 @param hist the histogram that contains the distribution 
	 @param n number of values requiered 
	 @param likelyhoodThreshold likelyhood threshold 
	 @param seriesProb probability of the obtained series 
	 @return  the random values 
	*/
	public final ArrayList<Double> GetRandomValues(Histogram hist, int n, double likelyhoodThreshold, DotNetHelpers.RefObject<Double> seriesProb)
	{
		ArrayList<Double> vals = new ArrayList<Double>();
		double val = -1.0;
		double prob;
		do
		{
			seriesProb.argValue = 0.0;
			vals.clear();
			for (int i = 0;i < n;i++)
			{
				prob = -1;
				while (prob < 0)
				{
					val = GetRandomValue(hist);
					prob = hist.GetFreq(val) / hist.TotFreqs;
				}
				if (seriesProb.argValue == 0)
				{
					seriesProb.argValue = prob;
				}
				else
				{
					seriesProb.argValue *= prob;
				}
				vals.add(val);
			}
		} while (seriesProb.argValue < likelyhoodThreshold);

		return vals;
	}

	/**  Get a secquence of n random values with the histogram distribution with a certain likelyhood 
	 @param hist the histogram that contains the distribution 
	 @param n number of values requiered 
	 @param likelyhoodThreshold likelyhood threshold 
	 @return  the random values 
	*/
	public final ArrayList<Double> GetRandomValues(Histogram hist, int n, double likelyhoodThreshold)
	{
		double seriesProb = 0;
		DotNetHelpers.RefObject<Double> tempRef_seriesProb = new DotNetHelpers.RefObject<Double>(seriesProb);
		List<Double> tempVar = GetRandomValues(hist, n, likelyhoodThreshold, tempRef_seriesProb);
		seriesProb = tempRef_seriesProb.argValue;
		return (ArrayList<Double>) tempVar;
	}

	/**  Get a set of grouped values (each group of nPerGroup elements) with their probabilities 
	 @param hist the histogram that contains the distribution 
	 @param nPerGroup number of values per group 
	 @param nGroupedValues number of grouped values 
	 @return  the set of grouped values 
	*/
	public final TreeMap<Integer, Double> GetRandomGroupedValues(Histogram hist, double nPerGroup, int nGroupedValues)
	{
		return GetRandomGroupedValues(hist, nPerGroup, nGroupedValues, -1.0);
	}

	/**  Get a set of grouped values (each group of nPerGroup elements) with their probabilities under a certain likekyhood threshold 
	 @param hist the histogram that contains the distribution 
	 @param nPerGroup number of values per group 
	 @param nGroupedValues number of grouped values 
	 @param likelyhoodThreshold likelyhood threshold 
	 @return  the set of grouped values 
	*/
	public final TreeMap<Integer, Double> GetRandomGroupedValues(Histogram hist, double nPerGroup, int nGroupedValues, double likelyhoodThreshold)
	{
		int nPerGroupInt = (int)nPerGroup;
		double diff = nPerGroup - (double)nPerGroupInt;

		TreeMap<Integer, Double> groupedVals = new TreeMap<Integer, Double>();
		double val;
		double prob = 0;
		int key;
		ArrayList<Double> seriesVals;
		for (int i = 0;i < nGroupedValues;i++)
		{
			DotNetHelpers.RefObject<Double> tempRef_prob = new DotNetHelpers.RefObject<Double>(prob);
			seriesVals = GetRandomValues(hist, nPerGroupInt + 1, likelyhoodThreshold, tempRef_prob);
		prob = tempRef_prob.argValue;
			val = func.Sum(seriesVals, 0, seriesVals.size() - 2);
			val += seriesVals.get(seriesVals.size() - 1) * diff;
			key = (int)Math.round(val);
			if (!groupedVals.containsKey(key))
			{
				groupedVals.put(key, prob);
			}
			else
			{
				groupedVals.put(key, groupedVals.get(key) + prob);
			}
		}
		return groupedVals;
	}

	private int BinarySearch(List<Double> values, double value)
	{
		return BinarySearch(values, value, 0, values.size() - 1);
	}

	private int BinarySearch(List<Double> values, double value, int low, int high)
	{
		if (high == -1 && (value <= values.get(low + 1)))
		{
			return low;
		}
		if (low == -1 && (high == 0 || value >= values.get(high - 1) && value >= values.get(high)))
		{
			return high;
		}
		if (low > high)
		{
			return -1;
		}
		if (high - low <= 1)
		{
			return (value < values.get(low))? low : high;
		}
		int i = low + (high - low) / 2;
		if ((value >= values.get(i - 1) && value <= values.get(i)))
		{
			return i;
		}
		if (value < values.get(i))
		{
			return BinarySearch(values, value, low, i - 1);
		}
		else
		{
			return BinarySearch(values, value, i + 1, high);
		}

	}

	//endregion

	//endregion

}