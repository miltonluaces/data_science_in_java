package NumCalc;

//region Imports

import ClassicalStat.*;
import java.util.*;

//endregion


/**  Class for density estimation by gaussian density functions composition 
*/
public class DensityEstim {

	//region Fields

	private Histogram histogram;
	private RootSearch rs;
	private NormalDist ns;
	private double s;
	private int maxIterations;
	private HashMap<Integer, Double> percentiles;
	private boolean mini;
	private boolean bisection;
	private ArrayList<Double> zerosDeriv2;
	private double minInt;
	private double maxInt;

	//endregion

	//region Constructor

	/**  Constructor 
	 @param s standard deviation for every gaussian function 
	 @param maxIterations maximum of iterations 
	 @param maxClasses maximum of classes 
	*/
	public DensityEstim(double s, int maxIterations, int maxClasses)
	{
		this.s = s;
		this.maxIterations = maxIterations;
		double epsilon = 0.005;
		double epsilonInterval = 0.05;
		histogram = new Histogram(maxClasses);
		percentiles = new HashMap<Integer, Double>();
		rs = new RootSearch(epsilon, epsilonInterval, maxIterations);
		ns = new NormalDist();
		mini = false;
		bisection = false;
		zerosDeriv2 = new ArrayList<Double>();
	}

	//endregion

	//region Properties, Setters and Getters

	/**  Maximum of iterations 
	*/
	public final int getMaxIterations()
	{
		return maxIterations;
	}
	public final void setMaxIterations(int value)
	{
		maxIterations = value;
	}

	/**  Histogram of the distribution 
	*/
	public final Histogram getHistogram()
	{
		return histogram;
	}
	public final void setHistogram(Histogram value)
	{
		histogram = value;
	}

	/**  Get the frequency of any value 
	 @param val the value 
	 @return  the frequency 
	*/
	public final double GetHistogramFreq(int val)
	{
		return histogram.GetFreq(val);
	}

	/**  Get percentile for any percentage 
	 @param perc the percentage 
	 @return  the percentile 
	*/
	public final double GetPercentile(double perc)
	{
		int key = GetKey(perc);
		return percentiles.get(key);
	}

	/**  Set a percentile for any value 
	 @param perc the percentage 
	 @param val the value 
	*/
	public final void SetPercentile(double perc, double val)
	{
		percentiles.put(GetKey(perc), val);
	}

	/**  if mini mode should be performed 
	*/
	public final boolean getMini()
	{
		return mini;
	}
	public final void setMini(boolean value)
	{
		mini = value;
	}

	/**  if bisection method should be performed (else montecarlo) 
	*/
	public final boolean getBisection()
	{
		return bisection;
	}
	public final void setBisection(boolean value)
	{
		bisection = value;
	}

	/**  Min value of the distribution 
	*/
	public final double getMin()
	{
		return histogram.Min;
	}

	/**  Max value of the distribution 
	*/
	public final double getMax()
	{
		return histogram.Max;
	}

	/**  Standard deviation for every gaussian function 
	*/
	public final double getS()
	{
		return s;
	}
	public final void setS(double value)
	{
		s = value;
	}

	//endregion

	//region Public Methods

	//region Load Methods

	/**  Load data from a list that represents a distribution (indexes = values, values = frequencies) 
	 @param meanFrec list with the distribution 
	*/
	public final void LoadDist(ArrayList<Double> meanFrec)
	{
		histogram.LoadDist(meanFrec);
	}

	/**  Load data from another histogram 
	 @param hist the other histogram 
	*/
	public final void LoadHist(Histogram hist)
	{
		histogram.LoadHist(hist);
	}

	public final void LoadDist(Map<Integer, Double> freqs)
	{
		histogram.LoadDist(freqs);
	}

	//endregion

	//region Probability Methods

	/**  Calculate probability for any value 
	 @param x the value 
	 @return  the calculated probability 
	*/
	public final double Probability(double x)
	{
		double p = 0.0;
		double totWeights = 0.0;
		double freq;
		for (double val : histogram.GetValues())
		{
			freq = histogram.GetFreq(val);
			p += freq * ns.pNorm(x, val, s);
			totWeights += freq;
		}
		if (totWeights == 0)
		{
			return 0.0;
		}
		return p / totWeights;
	}

	private double Quantile(double p)
	{
		return rs.SteffensenAccOneRoot(Probability, DerivProbability, Deriv2Probability, getMin(), getMax(), p);
	}

	/**  Calculate the percentile for any percentage 
	 @param p the percentage 
	 @return  the value for the percentile 
	*/
	public final double CalculatePercentile(double p)
	{
		if (histogram.size() == 1)
		{
			Iterator<Integer> enKeys = histogram.GetKeysEnumerator();
//C# TO JAVA CONVERTER TODO TASK: .NET iterators are only converted within the context of 'while' and 'for' loops:
			enKeys.MoveNext();
//C# TO JAVA CONVERTER TODO TASK: .NET iterators are only converted within the context of 'while' and 'for' loops:
			return (double)enKeys.Current;
		}

		double x = -1;
		int it = 1;
		double maxSteffensen = 0;
		//double maxSteffensen = 80;
		if (p >= maxSteffensen)
		{
			tangible.RefObject<Integer> tempRef_it = new tangible.RefObject<Integer>(it);
			x = rs.MonotoneBisection(Probability, true, 0.0, maxInt, p / 100.0, 0.005, tempRef_it, 120);
		it = tempRef_it.argValue;
		}
		else
		{
			try
			{
				x = rs.SteffensenAccOneRoot(Probability, DerivProbability, Deriv2Probability, 0.0, maxInt, p / 100.0, maxIterations);
			}
			catch (RuntimeException ex)
			{
				Trace.WriteLine("");
				x = -1;
			}

			if (x < 0 || Math.abs(Probability(x) - p / 100.0) > 0.01)
			{
				tangible.RefObject<Integer> tempRef_it2 = new tangible.RefObject<Integer>(it);
				x = rs.MonotoneBisection(Probability, true, 0.0, maxInt, p / 100.0, 0.005, tempRef_it2, 30);
			it = tempRef_it2.argValue;
			}
		}
		if (x < 0 || Math.abs(Probability(x) - p / 100.0) > 0.01)
		{
			x = -1;
		}
		//Console.WriteLine(x.ToString("0.00") + "\t" + p + "\t" + Probability(x).ToString("0.00%"));
		return x;
	}

	//endregion

	//region Set Percentiles

	/**  Set percentiles in variables 
	 @param allPercentiles if all percentiles must be calculated 
	*/
	public final int SetPercentiles(boolean allPercentiles)
	{
		percentiles.clear();
		if (mini)
		{
			SetPercentilesLess(allPercentiles);
			return 200;
		}
		else
		{
			maxInt = GetMaxInt();
			minInt = GetMinInt();
			if (allPercentiles)
			{
				for (int i = 800;i <= 1000;i++)
				{
					percentiles.put(i, 0);
				}
				int its = rs.SetAllValues(percentiles, Probability, 80, 100, minInt, maxInt, bisection);
				return its;
			}
			else
			{
				return 1;
			}
		}
	}


	/**  Set percentiles in case that data quantity is below the threshold 
	 @param allPercentiles if all percentiles must be calculated 
	*/
	private void SetPercentilesLess(boolean allPercentiles)
	{
		percentiles.clear();
		maxInt = GetMaxInt();
		minInt = GetMinInt();
		double x, xAnt, p;
		if (histogram.size() == 1)
		{
			Iterator<Integer> enKeys = histogram.GetKeysEnumerator();
//C# TO JAVA CONVERTER TODO TASK: .NET iterators are only converted within the context of 'while' and 'for' loops:
			enKeys.MoveNext();
			for (int i = 800;i < 1000;i++)
			{
				p = (double)i / 10.0;
				percentiles.put(GetKey(p), 0.0);
			}
			return;
		}
		xAnt = minInt;
		percentiles.put(GetKey(80.0), minInt);
		int it = 0;
		for (int i = 801;i < 999;i++)
		{
			p = (double)i / 10.0;
			try
			{
				tangible.RefObject<Integer> tempRef_it = new tangible.RefObject<Integer>(it);
				x = rs.MonotoneBisection(Probability, DerivProbability, xAnt, maxInt, p / 100.0, 0.0005, tempRef_it, 30);
			it = tempRef_it.argValue;
			}
			catch (java.lang.Exception e)
			{
				x = xAnt;
			}
			if (x < xAnt)
			{
				x = xAnt;
			}
			percentiles.put(GetKey(p), x);
			xAnt = x;
		}
		percentiles.put(GetKey(99.9), maxInt);
		percentiles.put(GetKey(100.0), maxInt);

	}

	/**  Set percentiles in case that data quantity is above the threshold 
	 @param allPercentiles if all percentiles must be calculated 
	*/
	private void SetPercentilesMore(boolean allPercentiles)
	{
		percentiles.clear();
		maxInt = GetMaxInt();
		minInt = GetMinInt();
		try
		{
			zerosDeriv2 = ZerosDeriv2();
		}
		catch (RuntimeException e)
		{
			throw new IllegalStateException(Strings.Cannot_calculate_second_derivative_zeros, e);
		}

		if (!allPercentiles)
		{
			return;
		}

		double xAnt = minInt;
		double x;
		int a = 0;
		int b = 1;
		double aInt, bInt;
		double p;
		percentiles.put(GetKey(80.0), minInt);
		for (int i = 801;i < 999;i++)
		{
			p = (double)i / 10.0;
			aInt = zerosDeriv2.get(a);
			bInt = zerosDeriv2.get(b);
			while (p / 100.0 > Probability(bInt) && b < zerosDeriv2.size() - 1)
			{
				a = b;
				b++;
				aInt = zerosDeriv2.get(a);
				bInt = zerosDeriv2.get(b);
			}
			try
			{
				x = rs.SteffensenAccOneRoot(Probability, DerivProbability, Deriv2Probability, aInt, bInt, p / 100.0, maxIterations);
			}
			catch (java.lang.Exception e)
			{
				//x = rs.MonotoneBisection(Probability, DerivProbability, aInt, bInt, p/100.0, 0.1, 100);
				x = xAnt;
			}
			aInt = x;
			if (x < xAnt)
			{
				x = xAnt;
			}
			if (x > maxInt)
			{
				x = maxInt;
			}
			percentiles.put(GetKey(p), x);
			xAnt = x;
		}
		percentiles.put(GetKey(99.9), maxInt);
		percentiles.put(GetKey(100.0), maxInt);
	}

	/**  Get minimum value of interval 
	 @return  the value 
	*/
	public final double GetMinInt()
	{
		int it = 0;
		tangible.RefObject<Integer> tempRef_it = new tangible.RefObject<Integer>(it);
		minInt = rs.MonotoneBisection(Probability, true, getMin(), maxInt, 0.8, 0.01, tempRef_it, 100);
	it = tempRef_it.argValue;
		return minInt;
	}

	/**  Get maximum value of interval 
	 @return  the value 
	*/
	public final double GetMaxInt()
	{
		int it = 0;
		double max = getMax() + getMax() * 0.1;
		while (Probability(max) < 0.99 && it < 10)
		{
			max += max * 0.1;
			it++;
		}
		return max;
	}

	/**  Set maximum interval (intermediate calculation) 
	*/
	public final void SetMaxInt()
	{
		maxInt = GetMaxInt();
	}

	//endregion

	//region Derivatives

	//region Weighted Normal Four Derivatives

	/** 
	 <p>La derivada primera de la función de probabilidad acumulada en x, es la función de densidad en x</p>
	 <p>Porque la función de probabilidad acumulada en x es la integral de la función de densidad hasta x</p>
	 <p>f' = Σ d(x)</p>
	 
	 @param x Valor en el cual se desea obtener la derivada de la función de probabilidad acumulada
	 @return  value of the derivative for x 
	*/
	public final double DerivProbability(double x)
	{
		double dens = 0.0;
		double freq;
		for (double val : histogram.GetValues())
		{
			freq = histogram.GetFreq(val);
			dens += freq * ns.pNormDeriv(x, (double)val, s);
		}
		return dens;
	}

	/** 
	 <p>Derivada segunda de la función de probabilidad acumulada en x.</p>
	 <p>f'' = Σ d'(x)  donde d'(x) = </p>
	 
	 @param x Valor en el cual se desea obtener la derivada segunda de la función de probabilidad acumulada
	 @return  value of the derivative for x 
	*/
	public final double Deriv2Probability(double x)
	{
		double dens = 0.0;
		double freq;
		for (int val : histogram.GetValues())
		{
			freq = histogram.GetFreq(val);
			dens += freq * ns.pNormDeriv2(x, (double)val, s);
		}
		return dens;
	}

	/** 
	 <p>Derivada tercera de la función de probabilidad acumulada en x.</p>
	 <p>f''' = Σ d''(x)  donde d'(x) = </p>
	 
	 @param x Valor en el cual se desea obtener la derivada tercera de la función de probabilidad acumulada
	 @return  value of the derivative for x 
	*/
	public final double Deriv3Probability(double x)
	{
		double dens = 0.0;
		double freq;
		for (int val : histogram.GetValues())
		{
			freq = histogram.GetFreq(val);
			dens += freq * ns.pNormDeriv3(x, (double)val, s);
		}
		return dens;
	}

	/** 
	 <p>Derivada cuarta de la función de probabilidad acumulada en x.</p>
	 <p>f'''' = Σ d'''(x)  donde d'(x) = </p>
	 
	 @param x Valor en el cual se desea obtener la derivada cuarta de la función de probabilidad acumulada
	 @return  value of the derivative for x 
	*/
	public final double Deriv4Probability(double x)
	{
		double dens = 0.0;
		double freq;
		for (int val : histogram.GetValues())
		{
			freq = histogram.GetFreq(val);
			dens += freq * ns.pNormDeriv4(x, (double)val, s);
		}
		return dens;
	}

	//endregion

	//region Zeros in Second Derivative

	//region General

	/**  Obtains a list of zeros in the second derivative of the function 
	 @return  List of zeros of second derivative 
	*/
	public final ArrayList<Double> ZerosDeriv2()
	{
		if (mini)
		{
			return ZerosDeriv2Less();
		}
		else
		{
			return ZerosDeriv2More();
		}
	}

	//endregion

	//region Less than 3 Frequencies

	private ArrayList<Double> ZerosDeriv2Less()
	{
		ArrayList<Double> zerosWDeriv2 = new ArrayList<Double>();
		zerosWDeriv2.add(minInt);
		ArrayList<Double> values = histogram.GetValues();
		double val1, val2, valInter;
		val1 = (double)values.get(0);
		val2 = (double)values.get(1);

		if (Deriv2Probability(val1) == 0 && val1 > minInt)
		{
			zerosWDeriv2.add(val1);
		}
		if (val2 == 0)
		{
			return zerosWDeriv2;
		}
		valInter = GetIntermediateZero(histogram.GetFreq(val1), histogram.GetFreq(val2), val1, val2);
		if (Deriv2Probability(valInter) == 0 && val1 > minInt)
		{
			zerosWDeriv2.add(valInter);
		}
		if (Deriv2Probability(val2) == 0 && val2 > minInt)
		{
			zerosWDeriv2.add(val2);
		}
		return zerosWDeriv2;
	}

	private double GetIntermediateZero(double w1, double w2, double m1, double m2)
	{
		return (2 * Math.log(w1 / w2) - Math.pow(m1, 2) + Math.pow(m2, 2)) / (2 * (m2 - m1));
	}

	//endregion

	//region More or equal 3 Frequencies

	private ArrayList<Double> ZerosDeriv2More()
	{
		zerosDeriv2.clear();
		zerosDeriv2.add(minInt);
		double interZero;
		double a, b, d2a, d2b;

		ArrayList<Double> values = histogram.GetValues();
		for (int i = 0;i < values.size() - 1;i++)
		{
			a = values.get(i);
			b = values.get(i + 1);

			d2a = Deriv2Probability(a);
			d2b = Deriv2Probability(b);

			if (d2a * d2b > -0.01)
			{
				a = b;
				continue;
			}
			interZero = rs.SteffensenAccOneRoot(Deriv2Probability, Deriv3Probability, Deriv4Probability, a, b, 0.0);
			if (interZero > minInt)
			{
				zerosDeriv2.add(interZero);
			}
		}
		zerosDeriv2.add(maxInt);
		return zerosDeriv2;
	}

	//endregion


	//endregion

	//endregion

	//endregion

	//region Private Methods

	private int GetKey(double val)
	{
		return (int)java.lang.Math.round(Math.round(val * Math.pow(10, 1)) / Math.pow(10, 1) * 10);
	}

	/**  Conversion key to value  
	 @param key the key 
	 @return  the corresponding value 
	*/
	public final double GetValue(int key)
	{
		return (double)key / 10.0;
	}

	private void AddPercentile(double perc, double val)
	{
		percentiles.put(GetKey(perc), val);
	}


	//endregion

}