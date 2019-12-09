package NumCalc;

import NumCalc.Properties.*;
import ClassicalStat.*;
import java.util.*;

//region Imports


//endregion

/**  Respository of miscellaneous functions 
*/
public class Functions
{

	//region Fields

	private Polynomic poly;
	private Statistics stat;
	private RandomGen rand;

	//endregion

	//region Constructor

	/**  Constructor  
	*/
	public Functions()
	{
		poly = new Polynomic();
		stat = new Statistics();
		rand = new RandomGen();
	}

	//endregion

	//region Distances

	/**  Euclidian distance 
	 @param Vector1 Vector of dobles that represents one point 
	 @param Vector2 Vector of dobles that represents the other point 
	 @return  calculated eucliedean distance 
	*/
	public final double DistEuclid(double[] Vector1, double[] Vector2)
	{
		if (Vector1.length != Vector2.length)
		{
			throw new RuntimeException(Strings.no_se_puede_calcular_distancia);
		}
		double dist = 0.0;
		for (int i = 0;i < Vector1.length;i++)
		{
			dist += Math.pow(Vector1[i] - Vector2[i], 2);
		}
		return Math.sqrt(dist);
	}

	//endregion

	//region Moving Average

	/**  Moving average algorithm 
	 @param values Lista de valores para los cuales se desea obtener la media movil
	 @param mw Numero de observaciones con las que se obtiene la media de cada punto
	 @param centered
		Si el numero de valores utilizados es mw/2 a la izquierda y otros tantos a la derecha. O si por el 
		contrario, son todos a la derecha del valor para el cual se calcula la media.
	 
	 @return  the time series modified 
	*/
	public final ArrayList<Double> MovingAverage(ArrayList<Double> values, int mw, boolean centered)
	{
		if (centered)
		{
			return MovingAverageCentered(values, mw);
		}
		else
		{
			return MovingAverage(values, mw);
		}
	}

	/**  Promedios moviles sin centrar
	 @param values Lista de valores para los cuales se desea obtener la media movil
	 @param mw Numero de observaciones a la derecha con las que se obtiene la media de cada punto
	 @return  the time series modified 
	*/
	private ArrayList<Double> MovingAverage(ArrayList<Double> values, int mw)
	{
		double average = 0.0;
		ArrayList<Double> movingAverage = new ArrayList<Double>();
		if (values.size() <= mw)
		{
			for (double val : values)
			{
				average += val;
			}
			for (double val : values)
			{
				movingAverage.add(average / values.size());
			}
		}
		else
		{
			//inicio
			for (int j = 0;j < mw;j++)
			{
				average += values.get(j);
			}
			for (int i = 0;i < mw;i++)
			{
				movingAverage.add(average / mw);
			}
			//medio
			for (int i = mw;i < values.size();i++)
			{
				average = 0.0;
				for (int j = i - mw + 1;j <= i;j++)
				{
					average += values.get(j);
				}
				movingAverage.add(average / mw);
			}
		}
		return movingAverage;
	}

	/**  Promedios moviles centrados
	 @param values Lista de valores para los cuales se desea obtener la media movil
	 @param mw Numero de observaciones centradas en el punto con las que se obtiene la media de cada punto
	 @return  the time series modified 
	*/
	private ArrayList<Double> MovingAverageCentered(ArrayList<Double> values, int mw)
	{
		ArrayList<Double> movingAverage = new ArrayList<Double>();
		int ini = (int)mw / 2;
		int fin;
		double average = 0.0;
		//inicio
		for (int j = 0;j < mw;j++)
		{
			average += values.get(j);
		}
		for (int i = 0;i < ini;i++)
		{
			movingAverage.add(average / mw);
		}
		//medio
		for (int i = ini;i < values.size() - ini;i++)
		{
			average = 0.0;
			fin = (mw % 2 == 0)? i + ini : i + ini + 1;
			for (int j = i - ini;j < fin;j++)
			{
				average += values.get(j);
			}
			movingAverage.add(average / mw);
		}
		//final
		average = 0.0;
		for (int j = values.size() - mw;j < values.size();j++)
		{
			average += values.get(j);
		}
		for (int i = values.size() - (mw - ini) - 1;i < values.size();i++)
		{
			movingAverage.add(average / mw);
			if (movingAverage.size() == values.size())
			{
				return movingAverage;
			}
		}
		return movingAverage;
	}

	//endregion

	//region Seasonal index

	/**  Calculate Seasonal index ( 
	 @param serie time series 
	 @param lag lag for autocorrelation 
	 @param cutoff cutoff for trimming 
	 @return  seasonal index 
	*/
	public final double CalcSeasonalIndex(ArrayList<Double> serie, int lag, double cutoff)
	{
		ArrayList<Double> ma = MovingAverage(serie, lag, true);
		ArrayList<Double> si = new ArrayList<Double>();
		for (int i = 0; i < serie.size(); i++)
		{
			si.add((serie.get(i) - ma.get(i)) / ma.get(i));
		}
		ArrayList<Double> trimSi = stat.Trim(si, cutoff, true);
		ArrayList<Double> ranges = new ArrayList<Double>();
		int ini = 0;
		int end = lag;
		while (end < si.size())
		{
			ranges.add(stat.Range(trimSi, ini, end));
			ini = end + 1;
			end = ini + 1 + lag;
		}
		return stat.Mean(ranges);
	}

	/**  Calculate Seasonal index ( 
	 @param serie time series 
	 @param lags list of lags for autocorrelation 
	 @param cutoff cutoff for trimming 
	 @return  list of seasonal indexes for each lag 
	*/
	public final ArrayList<Double> CalcSeasonalIndex(ArrayList<Double> serie, ArrayList<Integer> lags, double cutoff)
	{
		ArrayList<Double> seasIndexes = new ArrayList<Double>();
		double seasIndex;
		ArrayList<Double> ma, si, windSi, ranges;
		for (int lag : lags)
		{
			if (lag > serie.size() / 2.0)
			{
				continue;
			}
			ma = MovingAverage(serie, lag, true);
			si = new ArrayList<Double>();
			for (int i = 0; i < serie.size(); i++)
			{
				si.add((serie.get(i) - ma.get(i)));
			}
			windSi = stat.Windsor(si, cutoff, true);
			ranges = new ArrayList<Double>();
			int ini = 0;
			int end = lag;
			while (end < windSi.size())
			{
				ranges.add(stat.Range(windSi, ini, end));
				ini = end + 1;
				end = ini + 1 + lag;
			}

			seasIndex = stat.Mean(ranges) / stat.Mean(serie);
			seasIndexes.add(seasIndex);
		}
		return seasIndexes;
	}

	public final int CalcBestSeasonalLag(ArrayList<Double> serie, double cutoff, double minAutoCorr, double minAmplitude, int minLag, int maxLag)
	{
		ArrayList<Integer> lags = stat.CalculateEnoughAutocorrLags(serie, minAutoCorr, minLag, maxLag);
		if (lags == null || lags.isEmpty())
		{
			return -1;
		}
		if (lags.size() == 1)
		{
			return lags.get(0);
		}
		ArrayList<Double> sis = CalcSeasonalIndex(serie, lags, cutoff);
		int index = MaxIndex(sis);
		if (index == -1 || sis.get(index).compareTo(minAmplitude) < 0)
		{
			return -1;
		}
		else
		{
			return lags.get(index);
		}
	}

	//endregion

	//region Pendiente promedio

	/**  Average slope of a function, using angular coeff of quadratic means 
	 @param x independent variable 
	 @param y dependent variable (data) 
	 @return  the average slope 
	*/
	public final double AverageSlope(ArrayList<Double> x, ArrayList<Double> y)
	{
		ArrayList<Double> quadMins = poly.Regression(x, y, 1);
		return quadMins.get(1);
	}

	/**  Average slope of a function, using angular coeff of quadratic means 
	 @param vals data list 
	 @return  the average slope 
	*/
	public final double AverageSlope(ArrayList<Double> vals)
	{
		ArrayList<Double> xList = new ArrayList<Double>();
		for (int i = 0;i < vals.size();i++)
		{
			xList.add((double)i);
		}
		ArrayList<Double> quadMins = poly.Regression(xList, vals, 1);
		return quadMins.get(1);
	}

	//endregion

	//region Miscellaneous

	/**  Sum of values of a vector 
	 @param values the vector of values 
	 @return  the result sum 
	*/
	public final double Sum(List values)
	{
		return Sum(values, 0, values.size() - 1);
	}

	/**  Sum an interval of values of a vector 
	 @param values the vector of values 
	 @param start starting index 
	 @param end ending index 
	 @return  the result sum 
	*/
	public final double Sum(List values, int start, int end)
	{
		if (start < 0 || end > values.size())
		{
			throw new RuntimeException(Strings.Sum_out_of_interval);
		}
		double sum = 0.0;
		for (int i = start;i <= end;i++)
		{
			sum += (double)values.get(i);
		}
		return sum;
	}

	/**  Calculate min of an array of values 
	 @param valores the array of values 
	 @return  the min 
	*/
	public final double Min(double[] valores)
	{
		double min = 0.0;
		for (double val : valores)
		{
			if (val < min)
			{
				min = val;
			}
		}
		return min;
	}

	/**  Calculate max of an array of values 
	 @param valores the array of values 
	 @return  the max 
	*/
	public final double Max(ArrayList<Double> valores)
	{
		double max = 0.0;
		for (double val : valores)
		{
			if (val > max)
			{
				max = val;
			}
		}
		return max;
	}

	/**  Calculate max index of an array of values 
	 @param valores the array of values 
	 @return  the max index 
	*/
	public final int MaxIndex(ArrayList<Double> valores)
	{
		double max = 0.0;
		int index = -1;
		for (int i = 0;i < valores.size();i++)
		{
			if (valores.get(i).compareTo(max) > 0)
			{
				max = valores.get(i);
				index = i;
			}
		}
		return index;
	}

	/**  Normalize to positive data 
	 @param data data list 
	 @return  positive data list 
	*/
	public final ArrayList<Double> ConvertToPositive(ArrayList<Double> data)
	{
		double mindata = 0;
		ArrayList<Double> v = new ArrayList<Double>();

		for (int i = 0;i < data.size();i++)
		{
			if (data.get(i).compareTo(mindata) < 0)
			{
				mindata = data.get(i);
			}
		}
		mindata = Math.abs(mindata);
		for (int i = 0;i < data.size();i++)
		{
			v.add(data.get(i) + mindata);
		}
		return v;
	}

	/**  Calculate first intersection (from left to right) point between two functions 
	 @param arr1 first function data array 
	 @param arr2 second function data array 
	 @return  index of intersection 
	*/
	public final int GetFirstIntersectionIndex(double[] arr1, double[] arr2)
	{
		double lastSign = (int)Math.signum(arr1[0] - arr2[0]);
		double sign;
		for (int i = 1;i < arr1.length;i++)
		{
			sign = (int)Math.signum(arr1[i] - arr2[i]);
			if (sign != lastSign)
			{
				return i;
			}
		}
		if (lastSign > 0)
		{
			return -1;
		} //No se cortan, siendo el primer array siempre mayor
		else if (lastSign < 0)
		{
			return -2;
		} //No se cortan, siendo el segundo array siempre mayor
		return 0;
	}

	/**  Calculates series of increments or decrements (proportions) upon last value (multiplicative series) 
	 @param serie time series 
	 @return  multiplicative series 
	*/
	public final ArrayList<Double> GetMultiplicative(ArrayList<Double> serie)
	{
		ArrayList<Double> incFactorSerie = new ArrayList<Double>();
		incFactorSerie.add(1);
		for (int i = 1; i < serie.size(); i++)
		{
			incFactorSerie.add(serie.get(i) / serie.get(i - 1));
		}
		return incFactorSerie;
	}

	//"Solve" by Montecarlo newD2 (new data on data2) for a given newD1 (new datum on data1) and a r.
	public final double Solve(List<Double> data1, List<Double> data2, double r, double newD1, int iterations, int movingWindow, double maxInc)
	{
		rand.Reset();
		ArrayList<Double> data1Sim = new ArrayList<Double>(data1);
		ArrayList<Double> data2Sim = new ArrayList<Double>(data2);
		data1Sim.add(newD1);
		data2Sim.add(-1);
		double fcst = -1;
		double max = Max((ArrayList<Double>)data2);
		double d2Sim, rSim;
		double minDiff = Double.MAX_VALUE;
		for (int i = 0; i < iterations; i++)
		{
			d2Sim = (double)rand.NextInt(0, (int)(max * (1 + maxInc)));
			data2Sim.set(data2Sim.size() - 1, d2Sim);
			rSim = stat.R(data1Sim, data2Sim, data1Sim.size() - movingWindow - 1, data1Sim.size() - 1);
			double diff = Math.abs(r - rSim);
			if (diff < minDiff)
			{
				minDiff = diff;
				fcst = d2Sim;
			}
		}
		return fcst;
	}

	public final double Discretize(double value, double factor)
	{
		return Math.round((Math.round(value / factor) * factor) * Math.pow(10, 2)) / Math.pow(10, 2);
	}

	public final ArrayList<Double> Discretize(List<Double> values, double factor)
	{
		ArrayList<Double> discValues = new ArrayList<Double>();
		for (int i = 0; i < values.size(); i++)
		{
			discValues.add(Discretize(values.get(i), factor));
		}
		return discValues;
	}

	public final ArrayList<Double> SignFunction(List<Double> values)
	{
		ArrayList<Double> signFunction = new ArrayList<Double>();
		for (int i = 0; i < values.size(); i++)
		{
			if (values.get(i) < 0)
			{
				signFunction.add(-1);
			}
			else
			{
				signFunction.add(1);
			}
		}
		return signFunction;
	}

	//endregion

	//region Double lists operations

	/**  Substract the second list from the first list 
	 @param s1 first list 
	 @param s2 second list 
	 @param abs absolute value 
	 @param sort sort results 
	 @return  list with substraction result 
	*/
	public final ArrayList<Double> Substract(ArrayList<Double> s1, ArrayList<Double> s2, boolean abs, boolean sort)
	{
		if (s1.size() != s2.size())
		{
			throw new RuntimeException(Strings.s1_must_have_same_size_as_s2);
		}

		ArrayList<Double> res = new ArrayList<Double>();
		for (int i = 0;i < s1.size();i++)
		{
			if (abs)
			{
				res.add(Math.abs(s1.get(i)) - s2.get(i));
			}
			else
			{
				res.add(s1.get(i) - s2.get(i));
			}
		}
		if (sort)
		{
			Collections.sort(res);
		}
		return res;
	}

	/**  Calculate lower limit of the confidence interval (naive)
	 @param sortedList sorted list of data 
	 @param probability probability 
	 @param skipZeros if zeros shoud be skipped 
	 @return  lower limit of confidence interval 
	*/
	public final double GetInfLimitCI(ArrayList<Double> sortedList, double probability, boolean skipZeros)
	{
		if (probability < 0.0 || probability > 1.0)
		{
			throw new RuntimeException(Strings.probability_must_be_between_0_and_1);
		}
		int lastIndex = 0;
		int firstNoZero = 0;
		int nValues = 0;
		if (skipZeros)
		{
			for (int i = 0;i < sortedList.size();i++)
			{
				if (!sortedList.get(i).equals(0))
				{
					firstNoZero = i;
					nValues = sortedList.size() - firstNoZero;
					break;
				}
			}
		}
		else
		{
			nValues = sortedList.size() - 1;
		}
		lastIndex = firstNoZero + (int)Math.floor(nValues * probability); //+1
		if (lastIndex >= sortedList.size())
		{
			return Double.MAX_VALUE;
		}
		return sortedList.get(lastIndex);
	}

	/**  Filter outliers (naive) 
	 @param serie data serie 
	 @param model model data 
	 @param probability probability for the ci 
	 @param infLim lower limit 
	 @return  filtered data list 
	*/
	public final ArrayList<Double> FilterOutliers(ArrayList<Double> serie, ArrayList<Double> model, double probability, tangible.RefObject<Double> infLim)
	{
		ArrayList<Double> resid = Substract(serie, model, false, false);
		ArrayList<Double> sortedResid = new ArrayList<Double>();
		sortedResid.addAll(resid);
		Collections.sort(sortedResid);

		infLim.argValue = GetInfLimitCI(sortedResid, probability, true);
		//Console.WriteLine("InfLim = " + infLim);
		ArrayList<Double> filteredList = new ArrayList<Double>();
		double inf = Double.MAX_VALUE;
		for (int i = 0;i < resid.size();i++)
		{
			if (resid.get(i).compareTo(infLim.argValue) < 0)
			{
				filteredList.add(serie.get(i));
			}
			else
			{
				if (serie.get(i).compareTo(inf) < 0)
				{
					inf = serie.get(i);
				}
			}
			//Console.WriteLine("res = " + resid[i].ToString("0.0//") + " - serie = " + serie[i]);
		}
		infLim.argValue = inf;
		return filteredList;
	}

	//endregion

	//region Forgetting Function

	/**  Memory forgetting function for a serie 
	 @param totalPeriods total of series periods 
	 @param forgetInitPeriods start of forgetting period (from the present to the past) 
	 @param forgetEndPeriods end of forgetting period (from the present to the past) 
	 @param forgetEndProportion proportion of forgetting beyond the end (from the present to the past) 
	*/
	public final double[] CalcTempWeights(int totalPeriods, int forgetInitPeriods, int forgetEndPeriods, double forgetEndProportion)
	{
		double[] weights = new double[totalPeriods];
		for (int j = 0; j < totalPeriods; j++)
		{
			weights[j] = 1.00;
		}
		ApplyForgettingFunction(weights, forgetInitPeriods, forgetEndPeriods, forgetEndProportion);
		return weights;
	}


	public final double[] CalcTempWeights(int totalPeriods, int forgetInitPeriods, int forgetEndPeriods, double forgetEndProportion, List<Double> workingDays)
	{
		double[] weights = new double[workingDays.size()];
		for (int i = 0; i < workingDays.size(); i++)
		{
			weights[i] = workingDays.get(i);
		}
			ApplyForgettingFunction(weights, forgetInitPeriods, forgetEndPeriods, forgetEndProportion);
		return weights;
	}

	/**  Memory forgetting function  
	 @param weights array of weights to apply 
	 @param forgetInitPeriods start of forgetting period (from the present to the past) 
	 @param forgetEndPeriods en of forgetting period (from the present to the past) 
	 @param forgetEndProportion proportion of forgetting beyond the end (from the present to the past) 
	*/
	public final void ApplyForgettingFunction(double[] weights, int forgetInitPeriods, int forgetEndPeriods, double forgetEndProportion)
	{
		if (forgetInitPeriods == forgetEndPeriods)
		{
			return;
		}
		int forgetInitIndex = weights.length - forgetInitPeriods;
		int forgetEndIndex = weights.length - forgetEndPeriods;
		if (forgetInitIndex < 0)
		{
			forgetInitIndex = 0;
		}
		if (forgetEndIndex < 0)
		{
			forgetEndProportion = LinearFunction(0, forgetEndIndex, forgetEndProportion, forgetInitIndex, 1.0);
			forgetEndIndex = 0;
		}
		if (forgetInitIndex == forgetEndIndex)
		{
			return;
		}

		//constant final forgetting
		for (int i = 0;i < forgetEndIndex;i++)
		{
			weights[i] *= forgetEndProportion;
		}

		//variable ramp forgetting 
		double x1 = forgetEndIndex;
		double y1 = forgetEndProportion;
		double x2 = forgetInitIndex;
		double y2 = 1.0;

		double a = (y2 - y1) / (x2 - x1);
		double b = y1 - (x1 * (y2 - y1) / (x2 - x1));
		for (int x = forgetEndIndex;x < forgetInitIndex;x++)
		{
			weights[x] *= LinearFunction(x, x1, y1, x2, y2);
		}
	}

	/**  Linear function (using line equation) 
	 @param x independent variable 
	 @param x1 first defining point (x axis value) 
	 @param y1 first defining point (y axis value) 
	 @param x2 second defining point (x axis value) 
	 @param y2 second defining point (y axis value)
	 @return  functional value for x 
	*/
	public final double LinearFunction(double x, double x1, double y1, double x2, double y2)
	{
		double a = (y2 - y1) / (x2 - x1);
		double b = y1 - (x1 * (y2 - y1) / (x2 - x1));
		return a * x + b;
	}

	/**  Sigmoidal function between two points 
	 @param x independent variable 
	 @param x1 first defining point (x axis value) 
	 @param y1 first defining point (y axis value) 
	 @param x2 second defining point (x axis value) 
	 @param y2 second defining point (y axis value)
	 @return  functional value for x 
	*/
	public final double SigmoidalFunction(double x, double x1, double y1, double x2, double y2)
	{
		return 1.0 / (1.0 + Math.exp(-x));
	}

	//endregion

	//region Magnitude split


	/**  Discriminate magnitudes 
	 @param values values of the set or serie 
	 @param magnitudeBase base of logarithm for magnitude split 
	 @param normalized if values should be stored normalized or not 
	 @return  A list of lists of values where each list is stored in the index equal to exponent 
	*/
	public final ArrayList<ArrayList<Double>> MagnitudeSplit(ArrayList<Double> values, double magnitudeBase, boolean normalized)
	{
		ArrayList<ArrayList<Double>> magnitudes = new ArrayList<ArrayList<Double>>();
		for (int i = 0;i < 10;i++)
		{
			magnitudes.add(new ArrayList<Double>());
		}
		int exp;
		double denom;
		for (double value : values)
		{
			exp = (int)Math.log(value, magnitudeBase);
			if (normalized)
			{
				denom = Math.pow(magnitudeBase, exp);
				magnitudes.get(exp).add(Math.round(value / denom));
			}
			else
			{
				magnitudes.get(exp).add(value);
			}
		}
		while (magnitudes.get(magnitudes.size() - 1).isEmpty())
		{
			magnitudes.remove(magnitudes.size() - 1);
		}
		return magnitudes;
	}

	//endregion

	//region Random-Sequencies

	/**  Get distribution of sum of samples 
	 @param arrayToSample list to pick sample values (with repetition) 
	 @param sampleSize size of each sample 
	 @param nSamples number of samples 
	 @return  distribution of sample sums 
	*/
	public final HashMap<Double, Integer> GetSampleSumDistribution(List<Double> arrayToSample, int sampleSize, int nSamples)
	{
		HashMap<Double, Integer> sampSumDist = new HashMap<Double, Integer>();
		ArrayList<Double> sample;
		double sum;
		for (int i = 0; i < nSamples; i++)
		{
			sample = GetSample(arrayToSample, sampleSize);
			sum = Math.round((Sum(sample)) * Math.pow(10, 2)) / Math.pow(10, 2);
			if (!sampSumDist.containsKey(sum))
			{
				sampSumDist.put(sum, 0);
			}
			sampSumDist.get(sum)++;
		}
		return sampSumDist;
	}

	/**  Get list of sequences of sample sums different from a search value 
	 @param arrayToSample list to pick sample values (with repetition) 
	 @param sampleSize size of each sample 
	 @param nSamples number of samples 
	 @param searchValue the value to search 
	 @return  the list of sum sequences 
	*/
	public final ArrayList<Integer> GetNSequencesDifToK(List<Double> arrayToSample, int sampleSize, int nSamples, double searchValue)
	{
		searchValue = Math.round(searchValue * Math.pow(10, 2)) / Math.pow(10, 2);
		ArrayList<Integer> seqList = new ArrayList<Integer>();
		ArrayList<Double> sample;
		double sum;
		int seq = 0;
		for (int i = 0; i < nSamples; i++)
		{
			sample = GetSample(arrayToSample, sampleSize);
			sum = Math.round((Sum(sample)) * Math.pow(10, 2)) / Math.pow(10, 2);
			if (sum != searchValue)
			{
				seq++;
			}
			else
			{
				seqList.add(seq);
				seq = 0;
			}
		}
		if (seq != 0)
		{
			seqList.add(seq);
		}
		return seqList;
	}

	private ArrayList<Double> GetSample(List<Double> arrayToSample, int sampleSize)
	{
		int n = arrayToSample.size() - 1;
		ArrayList<Double> sample = new ArrayList<Double>();
		int index;
		for (int i = 0; i < sampleSize; i++)
		{
			index = rand.NextInt(0, n);
			sample.add(arrayToSample.get(index));
		}
		return sample;
	}

	//endregion

}