package ClassicalStat;

//region Imports

import java.util.*;


//endregion


/**  RandomProvider.  Provides random numbers of all data types in specified ranges.  It also contains a couple of methods from Normally (Gaussian) distributed random numbers and Exponentially distributed random numbers. 
*/
public class RandomGen
{

	//region Fields

	private Random rand;
	private double storedUniformDeviate;
	private boolean storedUniformDeviateIsGood = false;

	private static int seed = 0;
	private static boolean fast = true;

	//fast fields
	private static final double REAL_UNIT_INT = 1.0 / ((double)Integer.MAX_VALUE + 1.0);
	private static final double REAL_UNIT_UINT = 1.0 / ((double)Integer.MAX_VALUE + 1.0);
	private static final int Y = 842502087, Z = 357980759*10+1, W = 273326509;
	private int x, y, z, w;
	private int bitBuffer;
	private int bitMask = 1;


	//endregion

	//region Constructor

	/**  Constructor 
	*/
	public RandomGen()
	{
		Reset();
	}

	/**  Reset seed 
	*/
	public final void Reset()
	{
		if (fast)
		{
//C# TO JAVA CONVERTER WARNING: Unsigned integer types have no direct equivalent in Java:
//ORIGINAL LINE: x = (uint)(seed <= 0 ? Environment.TickCount : seed);
			x = (int)(seed <= 0 ? System.currentTimeMillis() : seed);
			y = Y;
			z = Z;
			w = W;
		}
		else
		{
			rand = new Random(seed <= 0 ? System.currentTimeMillis() : seed);
		}
	}

	/**  Force seed. 
	 @param newSeed seed for next random numbers set.
	*/
	public static void SetSeed(int newSeed)
	{
		seed = newSeed;
	}

	public static void SetFast(boolean isFast)
	{
		fast = isFast;
	}

	//endregion

	//region Distributions

	//region Uniform

//C# TO JAVA CONVERTER WARNING: Unsigned integer types have no direct equivalent in Java:
//ORIGINAL LINE: public uint NextUInt()
	public final int NextUInt()
	{
		if (fast)
		{
			return NextUIntFast();
		}
		synchronized (rand)
		{
//C# TO JAVA CONVERTER WARNING: Unsigned integer types have no direct equivalent in Java:
//ORIGINAL LINE: return (uint)Math.Abs(rand.Next());
			return (int)Math.abs(rand.nextInt());
		}
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
			throw new IllegalArgumentException("Strings.Max_must_be_greater_than_min");
		}
		if (fast)
		{
			return NextIntFast(min, max);
		}

		double nextDouble = 0.0;
		synchronized (rand)
		{
			nextDouble = rand.nextDouble();
		}
		int randInt = (int)java.lang.Math.round(min + nextDouble * (max - min));
		return randInt;
	}

	/**  Devuelve un double randómico en el rango [0,1) (cero inc) 
	 @return  result 
	*/
	public final double NextDouble()
	{
		if (fast)
		{
			return NextDoubleFast();
		}

		double nextDouble = 0.0;
		synchronized (rand)
		{
			nextDouble = rand.nextDouble();
		}
		return nextDouble;
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
			throw new IllegalArgumentException("Strings.Max_must_be_greater_than_min");
		}

		double nextDouble = 0.0;
		nextDouble = NextDouble();
		double randDbl = min + nextDouble * (max - min);
		return randDbl;
	}

	/**  Devuelve un booleano randómico 
	 @return  result 
	*/
	public final boolean NextBool()
	{
		if (fast)
		{
			return NextBoolFast();
		}
		double nextBoolean = 0.0;
		synchronized (rand)
		{
			nextBoolean = rand.nextInt(2);
		}
		return nextBoolean != 0;
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
		for (int i = 0; i < cantElementos; i++)
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

	/**  Devuelve sesgos randómicos positivos con media = 1, con distribucion Exponencial 
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

	/**  Devuelve sesgos randómicos positivos con media = m, con distribucion Exponencial 
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

		double val = NextDouble(0, hist.getTotFreqs());
		int valIndex = BinarySearch(hist.getAccumFreqValues(), val);
		int key = hist.GetKey(valIndex); //TODO: review
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
		for (int i = 0; i < n; i++)
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
			for (int i = 0; i < n; i++)
			{
				prob = -1;
				while (prob < 0)
				{
					val = GetRandomValue(hist);
					prob = hist.GetFreq(val) / hist.getTotFreqs();
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
	public final List<Double> GetRandomValues(Histogram hist, int n, double likelyhoodThreshold)
	{
		double seriesProb = 0;
		DotNetHelpers.RefObject<Double> tempRef_seriesProb = new DotNetHelpers.RefObject<Double>(seriesProb);
		List<Double> tempVar = GetRandomValues(hist, n, likelyhoodThreshold, tempRef_seriesProb);
		seriesProb = tempRef_seriesProb.argValue;
		return tempVar;
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
		for (int i = 0; i < nGroupedValues; i++)
		{
			DotNetHelpers.RefObject<Double> tempRef_prob = new DotNetHelpers.RefObject<Double>(prob);
			seriesVals = GetRandomValues(hist, nPerGroupInt + 1, likelyhoodThreshold, tempRef_prob);
		prob = tempRef_prob.argValue;
			val = 0.0;
			for (int j = 0; j <= seriesVals.size(); j++)
			{
				val += (double)seriesVals.get(j);
			}
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
			return (value < values.get(low)) ? low : high;
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

	//region Fast Methods

//C# TO JAVA CONVERTER WARNING: Unsigned integer types have no direct equivalent in Java:
//ORIGINAL LINE: private uint NextUIntFast()
	private int NextUIntFast()
	{
//C# TO JAVA CONVERTER WARNING: Unsigned integer types have no direct equivalent in Java:
//ORIGINAL LINE: uint t = (x ^ (x << 11));
		int t = (x ^ (x << 11));
		x = y;
		y = z;
		z = w;
//C# TO JAVA CONVERTER WARNING: The right shift operator was replaced by Java's logical right shift operator since the left operand was originally of an unsigned type, but you should confirm this replacement:
		return (w = (w ^ (w >>> 19)) ^ (t ^ (t >>> 8)));
	}

	private int NextIntFast()
	{
//C# TO JAVA CONVERTER WARNING: Unsigned integer types have no direct equivalent in Java:
//ORIGINAL LINE: uint t = (x ^ (x << 11));
		int t = (x ^ (x << 11));
		x = y;
		y = z;
		z = w;
//C# TO JAVA CONVERTER WARNING: The right shift operator was replaced by Java's logical right shift operator since the left operand was originally of an unsigned type, but you should confirm this replacement:
		return (int)(0x7FFFFFFF & (w = (w ^ (w >>> 19)) ^ (t ^ (t >>> 8))));
	}

	private int NextIntFast(int upper)
	{
		if (upper < 0)
		{
			throw new IllegalArgumentException("upperBound" + upper + "upperBound must be >=0");
		}
//C# TO JAVA CONVERTER WARNING: Unsigned integer types have no direct equivalent in Java:
//ORIGINAL LINE: uint t = (x ^ (x << 11));
		int t = (x ^ (x << 11));
		x = y;
		y = z;
		z = w;
//C# TO JAVA CONVERTER WARNING: The right shift operator was replaced by Java's logical right shift operator since the left operand was originally of an unsigned type, but you should confirm this replacement:
		return (int)((REAL_UNIT_INT * (int)(0x7FFFFFFF & (w = (w ^ (w >>> 19)) ^ (t ^ (t >>> 8))))) * upper);
	}

	private int NextIntFast(int lower, int upper)
	{
		if (lower > upper)
		{
			throw new IllegalArgumentException("upperBound" + upper + "upperBound must be >=lowerBound");
		}

//C# TO JAVA CONVERTER WARNING: Unsigned integer types have no direct equivalent in Java:
//ORIGINAL LINE: uint t = (x ^ (x << 11));
		int t = (x ^ (x << 11));
		x = y;
		y = z;
		z = w;
		int range = upper - lower;
		if (range < 0)
		{
//C# TO JAVA CONVERTER WARNING: The right shift operator was replaced by Java's logical right shift operator since the left operand was originally of an unsigned type, but you should confirm this replacement:
			return lower + (int)((REAL_UNIT_UINT * (double)(w = (w ^ (w >>> 19)) ^ (t ^ (t >>> 8)))) * (double)((long)upper - (long)lower));
		}
//C# TO JAVA CONVERTER WARNING: The right shift operator was replaced by Java's logical right shift operator since the left operand was originally of an unsigned type, but you should confirm this replacement:
		return lower + (int)((REAL_UNIT_INT * (double)(int)(0x7FFFFFFF & (w = (w ^ (w >>> 19)) ^ (t ^ (t >>> 8))))) * (double)range);
	}

	private double NextDoubleFast()
	{
//C# TO JAVA CONVERTER WARNING: Unsigned integer types have no direct equivalent in Java:
//ORIGINAL LINE: uint t = (x ^ (x << 11));
		int t = (x ^ (x << 11));
		x = y;
		y = z;
		z = w;
//C# TO JAVA CONVERTER WARNING: The right shift operator was replaced by Java's logical right shift operator since the left operand was originally of an unsigned type, but you should confirm this replacement:
		return (REAL_UNIT_INT * (int)(0x7FFFFFFF & (w = (w ^ (w >>> 19)) ^ (t ^ (t >>> 8)))));
	}

	private boolean NextBoolFast()
	{
		if (bitMask == 1)
		{
//C# TO JAVA CONVERTER WARNING: Unsigned integer types have no direct equivalent in Java:
//ORIGINAL LINE: uint t = (x ^ (x << 11));
			int t = (x ^ (x << 11));
			x = y;
			y = z;
			z = w;
//C# TO JAVA CONVERTER WARNING: The right shift operator was replaced by Java's logical right shift operator since the left operand was originally of an unsigned type, but you should confirm this replacement:
			bitBuffer = w = (w ^ (w >>> 19)) ^ (t ^ (t >>> 8));

			bitMask = (int)0x80000000;
			return (bitBuffer & bitMask) == 0;
		}

//C# TO JAVA CONVERTER WARNING: The right shift operator was replaced by Java's logical right shift operator since the left operand was originally of an unsigned type, but you should confirm this replacement:
		return (bitBuffer & (bitMask >>>= 1)) == 0;
	}

	//not used yet
//C# TO JAVA CONVERTER WARNING: Unsigned integer types have no direct equivalent in Java:
//ORIGINAL LINE: private void NextBytesFast(byte[] buffer)
	private void NextBytesFast(byte[] buffer)
	{

//C# TO JAVA CONVERTER WARNING: Unsigned integer types have no direct equivalent in Java:
//ORIGINAL LINE: uint x = this.x, y = this.y, z = this.z, w = this.w;
		int x = this.x, y = this.y, z = this.z, w = this.w;
		int i = 0;
//C# TO JAVA CONVERTER WARNING: Unsigned integer types have no direct equivalent in Java:
//ORIGINAL LINE: uint t;
		int t;
		for (int bound = buffer.length - 3; i < bound;)
		{
			t = (x ^ (x << 11));
			x = y;
			y = z;
			z = w;
//C# TO JAVA CONVERTER WARNING: The right shift operator was replaced by Java's logical right shift operator since the left operand was originally of an unsigned type, but you should confirm this replacement:
			w = (w ^ (w >>> 19)) ^ (t ^ (t >>> 8));

//C# TO JAVA CONVERTER WARNING: Unsigned integer types have no direct equivalent in Java:
//ORIGINAL LINE: buffer[i++] = (byte)w;
			buffer[i++] = (byte)w;
//C# TO JAVA CONVERTER WARNING: The right shift operator was replaced by Java's logical right shift operator since the left operand was originally of an unsigned type, but you should confirm this replacement:
//C# TO JAVA CONVERTER WARNING: Unsigned integer types have no direct equivalent in Java:
//ORIGINAL LINE: buffer[i++] = (byte)(w >> 8);
			buffer[i++] = (byte)(w >>> 8);
//C# TO JAVA CONVERTER WARNING: The right shift operator was replaced by Java's logical right shift operator since the left operand was originally of an unsigned type, but you should confirm this replacement:
//C# TO JAVA CONVERTER WARNING: Unsigned integer types have no direct equivalent in Java:
//ORIGINAL LINE: buffer[i++] = (byte)(w >> 16);
			buffer[i++] = (byte)(w >>> 16);
//C# TO JAVA CONVERTER WARNING: The right shift operator was replaced by Java's logical right shift operator since the left operand was originally of an unsigned type, but you should confirm this replacement:
//C# TO JAVA CONVERTER WARNING: Unsigned integer types have no direct equivalent in Java:
//ORIGINAL LINE: buffer[i++] = (byte)(w >> 24);
			buffer[i++] = (byte)(w >>> 24);
		}

		if (i < buffer.length)
		{
			t = (x ^ (x << 11));
			x = y;
			y = z;
			z = w;
//C# TO JAVA CONVERTER WARNING: The right shift operator was replaced by Java's logical right shift operator since the left operand was originally of an unsigned type, but you should confirm this replacement:
			w = (w ^ (w >>> 19)) ^ (t ^ (t >>> 8));

//C# TO JAVA CONVERTER WARNING: Unsigned integer types have no direct equivalent in Java:
//ORIGINAL LINE: buffer[i++] = (byte)w;
			buffer[i++] = (byte)w;
			if (i < buffer.length)
			{
//C# TO JAVA CONVERTER WARNING: The right shift operator was replaced by Java's logical right shift operator since the left operand was originally of an unsigned type, but you should confirm this replacement:
//C# TO JAVA CONVERTER WARNING: Unsigned integer types have no direct equivalent in Java:
//ORIGINAL LINE: buffer[i++] = (byte)(w >> 8);
				buffer[i++] = (byte)(w >>> 8);
				if (i < buffer.length)
				{
//C# TO JAVA CONVERTER WARNING: The right shift operator was replaced by Java's logical right shift operator since the left operand was originally of an unsigned type, but you should confirm this replacement:
//C# TO JAVA CONVERTER WARNING: Unsigned integer types have no direct equivalent in Java:
//ORIGINAL LINE: buffer[i++] = (byte)(w >> 16);
					buffer[i++] = (byte)(w >>> 16);
					if (i < buffer.length)
					{
//C# TO JAVA CONVERTER WARNING: The right shift operator was replaced by Java's logical right shift operator since the left operand was originally of an unsigned type, but you should confirm this replacement:
//C# TO JAVA CONVERTER WARNING: Unsigned integer types have no direct equivalent in Java:
//ORIGINAL LINE: buffer[i] = (byte)(w >> 24);
						buffer[i] = (byte)(w >>> 24);
					}
				}
			}
		}
		this.x = x;
		this.y = y;
		this.z = z;
		this.w = w;
	}

	//endregion

}