package ClassicalStat;

import java.util.*;

/**  Container class for statistical distributions 
*/
public class Histogram
{

	//region Fields

	private ECCSortedDictionary<Integer, Double> freqs;
	private ECCSortedDictionary<Integer, Double> accumFreqs;
	private ArrayList<Integer> accumFreqKeys;
	private ArrayList<Double> accumFreqValues;
	private double totalFreqs;
	private int maxClasses;
	private double min;
	private double max;
	private double mean;
	private double stDev;
	private boolean normalized;
	private int loads;
	private RandomGen rand;
	private double width;
	private int maxKey;

	private double total;
	private double error;
	private boolean statsCalculated;

	private HashMap<Double, Double> cacheProbabilities;

	//endregion

	//region Constructors

	/**  Constructor 
	*/
	public Histogram()
	{
		this.maxClasses = Integer.MAX_VALUE;
		freqs = new ECCSortedDictionary<Integer, Double>();
		accumFreqs = new ECCSortedDictionary<Integer, Double>();
		accumFreqKeys = new ArrayList<Integer>();
		accumFreqValues = new ArrayList<Double>();
		cacheProbabilities = new HashMap<Double, Double>();
		totalFreqs = -1.0;
		width = -1;
		mean = -1;
		Clear();
		rand = new RandomGen();
		rand.Reset();
		total = 0;
		error = 0;
		statsCalculated = false;
	}

	/**  Constructor 
	 @param maxClasses maximum of hystogram 
	*/
	public Histogram(int maxClasses)
	{
		this();
		this.maxClasses = maxClasses;
	}


	//endregion

	//region Properties

	/**  maximum number of classes allowed 
	*/
	public final int getMaxClasses()
	{
		return maxClasses;
	}

	/**  min value 
	*/
	public final double getMin()
	{
		return min;
	}

	/**  max value 
	*/
	public final void setMax(double value)
	{
		max = value;
	}
	public final double getMax()
	{
		return max;
	}

	/**  values range 
	*/
	public final double getRange()
	{
		return max - min;
	}

	/**  column width 
	*/
	public final double getWidth()
	{
		if (width == -1)
		{
			width = (getRange() + 1) / maxKey;
		}
		return width;
	}

	/**  if values are normalized or not 
	*/
	public final boolean getNormalized()
	{
		return normalized;
	}

	/**  quantity of values 
	*/
	public final int getCount()
	{
		return freqs.size();
	}

	/**  number of times when some data has been loaded 
	*/
	public final int getLoads()
	{
		return loads;
	}

	/**  sorted collection of frequencies 
	*/
	public final Map<Integer, Double> getFreqs()
	{
		return freqs;
	}
	public final void setFreqs(Map<Integer, Double> value)
	{
		freqs = (ECCSortedDictionary<Integer, Double>) value;
	}

	/**  sorted collection of accumulated frequencies 
	*/
	public final ECCSortedDictionary<Integer, Double> getAccumFreqs()
	{
		return accumFreqs;
	}

	/**  list of accumulated frequencies keys 
	*/
	public final ArrayList<Integer> getAccumFreqKeys()
	{
		return accumFreqKeys;
	}

	/**  list of accumulated frequencies values 
	*/
	public final ArrayList<Double> getAccumFreqValues()
	{
		return accumFreqValues;
	}

	/**  total of all frequencies 
	*/
	public final double getTotFreqs()
	{
		if (totalFreqs == -1)
		{
			CalculateStatistics();
		}
		return totalFreqs;
	}

	/**  Distribution mean 
	*/
	public final double getMean()
	{
		if (mean == -1)
		{
			CalculateStatistics();
		}
		return mean;
	}

	/**  Distribution standard deviation 
	*/
	public final double getStDev()
	{
		if (stDev == -1)
		{
			CalculateStatistics();
		}
		return stDev;
	}


	/**  Relative error 
	*/
	public final double getRelError()
	{
		return error / total;
	}

	/**  Non zero freqs from the histogram 
	*/
	public final int getNonZeroFreqs()
	{
		return freqs.size();
	}

	//endregion

	//region Public Methods

	//region Getters

	/**  get the frequency of a certain value 
	 @param value the value 
	 @return  the frequency 
	*/
	public final double GetFreq(double value)
	{
		int key = GetKey(value);
		if (!freqs.containsKey(key))
		{
			return -1;
		}
		return freqs.get(key);
	}

	/**  get the accumulated frequency of a certain value 
	 @param value the value 
	 @return  the accumulated frequency 
	*/
	public final double GetAccumFreq(double value)
	{
		int key = GetKey(value);
		if (!accumFreqs.containsKey(key))
		{
			return -1;
		}
		return accumFreqs.get(key);
	}

	/**  If the histogram contains a certain value 
	 @param value the value 
	 @return  if exists 
	*/
	public final boolean Contains(double value)
	{
		int key = GetKey(value);
		return freqs.containsKey(key);
	}

	/**  Get normalized value 
	 @param value the value 
	 @return  the normalized value 
	*/
	public final double GetNormValue(double value)
	{
		if (!normalized)
		{
			return value;
		}
		return GetValue(GetKey(value));
	}

	//endregion

	//region Load Data

	/**  Load a collection of values on the hystogram 
	 @param values the collection of values 
	*/
	public final void LoadData(ArrayList<Double> values)
	{
		SetDataMinMax(values);
		normalized = (getRange() > maxClasses);
		for (double value : values)
		{
			AddValue(value);
		}
		loads++;
	}

	/**  Load a collection that represent frequencies (indexes represent values) 
	 @param freqs the collection of frequencies 
	*/
	public final void LoadDist(ArrayList<Double> freqs)
	{
		SetDistMinMax(freqs);
		normalized = (getRange() > maxClasses);
		for (int i = (int)min;i <= (int)max;i++)
		{
			AddValue((double)i, freqs.get(i));
		}
		loads++;
	}

	/**  Load a sorted dictionary of frequencies 
	 @param sortFreqs the dictionary 
	*/
	public final void LoadDist(ECCSortedDictionary<Integer, Double> sortFreqs)
	{
		SetDistMinMax(sortFreqs);
		normalized = (getRange() > maxClasses);
		if (!normalized)
		{
			freqs = sortFreqs;
		}
		else
		{
			for (int key : sortFreqs.keySet())
			{
				AddValue((double)key, sortFreqs.get(key));
			}
		}
		loads++;
	}

	/**  Load a sorted dictionary with min and max already calculated 
	 @param sortFreqs the dictionary 
	 @param min calculated min 
	 @param max calculated max 
	*/
	public final void LoadDist(Map<Integer, Double> sortFreqs, double min, double max)
	{
		if (min < this.min)
		{
			this.min = min;
		}
		if (max > this.max)
		{
			this.max = max;
		}
		normalized = (getRange() > maxClasses);
		for (int key : sortFreqs.keySet())
		{
				AddValue((double)key, sortFreqs.get(key)); //TODO: revisar si debe ser key o value el primer parametro
		}
		loads++;
	}

	public final void LoadDist(Map<Integer, Double> sortFreqs)
	{
		SetDistMinMax(sortFreqs);
		normalized = (getRange() > maxClasses);
		for (int key : sortFreqs.keySet())
		{
			AddValue((double)key, sortFreqs.get(key)); //TODO: revisar si debe ser key o value el primer parametro
		}
		loads++;
	}

	/**  Load an histogram 
	 @param hist the histogram 
	*/
	public final void LoadHist(Histogram hist)
	{
		this.min = hist.getMin();
		this.max = hist.getMax();
		this.totalFreqs = hist.getTotFreqs();
		normalized = (getRange() > maxClasses);

		ArrayList<Integer> sortedkeys = new ArrayList<Integer>(hist.freqs.keySet());
		Collections.sort(sortedkeys);
		for (int key : sortedkeys)
		{
			AddValue(hist.GetValue(key), hist.freqs.get(key));
		}
	}

	/**  Load two convoluted histograms 
	 @param h1 the histogram 1
	 @param h2 the histogram 2
	*/
	public final void LoadHist(Histogram h1, Histogram h2)
	{
		ArrayList<Double> freqs = h1.GetFreqs();
		freqs.addAll(h2.GetFreqs());
		this.LoadDist(freqs);
	}

	public final void LoadData(ArrayList<Double> values, List<Double> weights)
	{
		SetDataMinMax(values);
		normalized = (getRange() > maxClasses);
		for (int i = 0; i < values.size();i++)
		{
			AddValue(values.get(i), weights.get(i));
		}
		loads++;
	}


	/**  Calculate statistics 
	*/
	public final void CalculateStatistics()
	{
		if (statsCalculated)
		{
			return;
		}

		double accum = 0.0;
		double sum = 0.0;
		double sumSq = 0.0;
		double val, freq;

		ArrayList<Integer> sortedkeys = new ArrayList<Integer>(freqs.keySet());
		Collections.sort(sortedkeys);
		for (int key : sortedkeys)
		{
			val = GetValue(key);
			freq = freqs.get(key);
			accum += freq;
			accumFreqKeys.add(key);
			accumFreqValues.add(accum);
			accumFreqs.Add(key, accum);
			sum += (val * freq);
			sumSq += (val * val * freq);
		}
		totalFreqs = accum;
		mean = sum / (double)totalFreqs;
		double sqSum = sum * sum;
		if (totalFreqs == 0 || freqs.size() < 2)
		{
			stDev = 0;
		}
		else
		{
			stDev = Math.sqrt((sumSq - (sqSum / totalFreqs)) / (totalFreqs));
		}
		if (Double.isNaN(stDev))
		{
			throw new RuntimeException("Strings.StDev_cannot_be_negative");
		}
		statsCalculated = true;
	}

	/**  Clear the hystogram 
	*/
	public final void Clear()
	{
		freqs.clear();
		min = Double.MAX_VALUE;
		max = -Double.MAX_VALUE;
		normalized = false;
	}

	//endregion

	//region To Lists

	/**  Get a list of hystogram values 
	 @return  a list of doubles 
	*/
	public final ArrayList<Double> GetValues()
	{
		ArrayList<Double> values = new ArrayList<Double>();
		for (int key : freqs.keySet())
		{
			values.add(GetValue(key));
		}
		return values;
	}

	/**  Get a list of hystogram frequencies 
	 @return  a list of values 
	*/
	public final ArrayList<Double> GetFreqs()
	{
		ArrayList<Double> frequencies = new ArrayList<Double>();
		for (int key : freqs.keySet())
		{
			frequencies.add(freqs.get(key));
		}
		return frequencies;
	}

	/**  Get an array of normalized frequencies 
	 @return  an array of doubles 
	*/
	public final double[] GetArrayNormFreqs()
	{
		double[] arrayFreqs = new double[(int)max + 1];
		for (int key : freqs.keySet())
		{
			arrayFreqs[key] = freqs.get(key);
		}
		return arrayFreqs;
	}

	/**  Get an enumerator of values and frequencies 
	 @return  the enumerator 
	*/
	public final Iterator<Integer> GetKeysEnumerator()
	{
		return freqs.keySet().iterator();
	}

	/**  Get list of keys 
	 @return  list of keys 
	*/
	public final ArrayList<Integer> GetKeys()
	{
		ArrayList<Integer> keys = new ArrayList<Integer>();
		for (int key : freqs.keySet())
		{
			keys.add(key);
		}
		return keys;
	}

	/**  Get a list of frequencies 
	 @return  a list of doubles 
	*/
	public final ArrayList<Double> GetFreqsList()
	{
		ArrayList<Double> freqsList = new ArrayList<Double>();
		for (int i = 0;i < maxClasses;i++)
		{
			if (freqs.containsKey(i))
			{
				freqsList.add(freqs.get(i));
			}
			else
			{
				freqsList.add(0.0);
			}
		}
		return freqsList;
	}

	/**  Discrete probability (classical definition) 
	 @param value value to calculate probability 
	 @return  probability value 
	*/
	public final double ProbabilityByValue(double value)
	{
		//if(cacheProbabilities.ContainsKey(value)) { return cacheProbabilities[value]; } 
		double freqValue = GetFreq(value);
		if (freqValue <= 0)
		{
			return 0;
		}
		double prob = freqValue / getTotFreqs();
		//cacheProbabilities.Add(value, prob);
		return prob;

	}

	/**  Discrete probability (classical definition) 
	 @param key key to calculate probability 
	 @return  probability value 
	*/
	public final double ProbabilityByKey(int key)
	{
		if (cacheProbabilities.containsKey(key))
		{
			return cacheProbabilities.get(key);
		}
		double freqValue = freqs.get(key);
		if (freqValue <= 0)
		{
			return 0;
		}
		double prob = freqValue / getTotFreqs();
		cacheProbabilities.putIfAbsent((double) key, prob);
		return prob;
	}

	/**  Get a list of raw data from the histogram 
	 @return  the raw data list 
	*/
	public final ArrayList<Double> GetRawData()
	{
		double normFactor = 1000 / getTotFreqs();
		;
		ArrayList<Double> values = GetValues();
		ArrayList<Double> freqs = GetFreqs();
		ArrayList<Double> rawData = new ArrayList<Double>();
		for (int i = 0; i < values.size(); i++)
		{
			for (int j = 0; j < (int)(Math.round(freqs.get(i) * normFactor)); j++)
			{
				rawData.add(values.get(i));
			}
		}
		return rawData;
	}

	//endregion

	//endregion

	//region Private Methods

	/**  get the int normalized value of a value 
	 @param value double value 
	 @param total total value 
	 @param error error value 
	 @return  the key 
	*/
	public final int GetKey(double value, DotNetHelpers.RefObject<Double> total, DotNetHelpers.RefObject<Double> error)
	{
		if (normalized)
		{
			value = Normalize(value);
		}
		int round = (int)Math.round(value);
		total.argValue = Math.max(value, round);
		error.argValue = Math.abs(value - round);
		return round;
	}

	/**  get the int normalized value of a value 
	 @param value double value 
	 @return  the key 
	*/
	public final int GetKey(double value)
	{
		if (normalized)
		{
			value = Normalize(value);
		}
		int round = (int)Math.round(value);
		return round;
	}

	/**  get the double unnormalized value of a key 
	 @param key int key value 
	 @return  the double value 
	*/
	public final double GetValue(int key)
	{
		if (normalized)
		{
			return UnNormalize((double)key);
		}
		return (double)key;
	}

	/**  get the int normalized value of a value for convolutions
	 @param value double value 
	 @param n number of convolutions 
	 @return  the key 
	*/
	public final int GetConvKey(double value, int n)
	{
		if (normalized)
		{
			value = NormalizeConv(value, n);
		}
		int round = (int)Math.round(value);
		return round;
	}

	/**  get the double unnormalized value of a key for convolutions
	 @param key int key value 
	 @param n number of convolutions 
	 @return  the double value 
	*/
	public final double GetConvValue(int key, int n)
	{
		if (normalized)
		{
			return UnNormalizeConv((double)key, n);
		}
		return (double)key;
	}

	private void AddValue(double value)
	{
		double tot = -1;
		double err = -1;
		DotNetHelpers.RefObject<Double> tempRef_tot = new DotNetHelpers.RefObject<Double>(tot);
		DotNetHelpers.RefObject<Double> tempRef_err = new DotNetHelpers.RefObject<Double>(err);
		int key = GetKey(value, tempRef_tot, tempRef_err);
	err = tempRef_err.argValue;
	tot = tempRef_tot.argValue;
		total += tot;
		error += err;
		if (freqs.containsKey(key))
		{
			freqs.setItem(key, freqs.get(key) + 1.0);
		}
		else
		{
			freqs.Add(key, 1.0);
			if (key > maxKey)
			{
				maxKey = key;
			}
		}
	}

	private void AddValue(double value, double freq)
	{
		if (freq == 0)
		{
			return;
		}
		double tot = -1;
		double err = -1;
		DotNetHelpers.RefObject<Double> tempRef_tot = new DotNetHelpers.RefObject<Double>(tot);
		DotNetHelpers.RefObject<Double> tempRef_err = new DotNetHelpers.RefObject<Double>(err);
		int key = GetKey(value, tempRef_tot, tempRef_err);
	err = tempRef_err.argValue;
	tot = tempRef_tot.argValue;
		total += tot * freq;
		error += err * freq;
		if (freqs.containsKey(key))
		{
			freqs.setItem(key, freqs.get(key) + freq);
		}
		else
		{
			freqs.Add(key, freq);
			if (key > maxKey)
			{
				maxKey = key;
			}
		}
	}

	private void SetDataMinMax(ArrayList<Double> values)
	{
		for (int i = 0; i < values.size();i++)
		{
			if (values.get(i).compareTo(0.0) < 0)
			{
				values.set(i, 0.0);
			}
			if (values.get(i).compareTo(min) < 0)
			{
				min = values.get(i);
			}
			if (values.get(i).compareTo(max) > 0)
			{
				max = values.get(i);
			}
		}
	}

	private void SetDistMinMax(ArrayList<Double> freqs)
	{
		for (int i = 0;i < freqs.size();i++)
		{
			if (freqs.get(i).equals(0))
			{
				continue;
			}
			if (min == Double.MAX_VALUE)
			{
				min = (double)i;
			}
			max = (double)i;
		}
	}

	private void SetDistMinMax(Map<Integer, Double> freqs)
	{
		for (int key : freqs.keySet())
		{
			if (min == Double.MAX_VALUE)
			{
				min = (double)key;
			}
			max = (double)key;
		}
	}

	private double Normalize(double value)
	{
		return ((value - min) / getRange()) * maxClasses;
	}

	private double UnNormalize(double value)
	{
		return (value * getRange()) / maxClasses + min;
	}

	private double NormalizeConv(double value, int n)
	{
		return ((value - min * n) / getRange()) * maxClasses;
	}

	private double UnNormalizeConv(double value, int n)
	{
		return (value * getRange()) / maxClasses + min * n;
	}

	//endregion
}