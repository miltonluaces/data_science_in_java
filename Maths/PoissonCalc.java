package NumCalc;

//region Imports

import ClassicalStat.*;
import java.util.*;

//endregion


/**  Class for any Poisson distribution calculations 
*/
public class PoissonCalc
{

	//region Fields

	private ArrayList<Double> data;
	private double loadFactor;

	private ArrayList<Double>[] dataSplit;
	private double[] Lambda;
	private double[] Alpha;
	private double alpha;
	private double lambda;

	private DiscreteDist dd;
	private Statistics stat;
	private RootSearch rs;
	private RandomGen rg;
	private Functions func;

	private double max;
	private double epsilon;
	private int maxIt;

	private MethodType method = MethodType.values()[0];

	private double[] weights;
	private int forgetInitPeriods;
	private int forgetEndPeriods;
	private double forgetEndProportion;
	private double firstNonZero;


	//endregion

	//region Constructor

	/**  Constructor 
	 @param epsilon minimum value for approximation 
	 @param maxIt maximum number of iterations 
	 @param forgetInitPeriods start of forgetting period (from the present to the past) 
	 @param forgetEndPeriods en of forgetting period (from the present to the past) 
	 @param forgetEndProportion proportion of forgetting beyond the end (from the present to the past) 
	*/
	public PoissonCalc(double epsilon, int maxIt, int forgetInitPeriods, int forgetEndPeriods, double forgetEndProportion)
	{
		this.dd = new DiscreteDist();
		this.stat = new Statistics();
		this.epsilon = epsilon;
		this.maxIt = maxIt;
		this.rs = new RootSearch(epsilon, -1, maxIt);
		this.rg = new RandomGen();
		this.func = new Functions();
		this.method = MethodType.Poisson;
		this.forgetInitPeriods = forgetInitPeriods;
		this.forgetEndPeriods = forgetEndPeriods;
		this.forgetEndProportion = forgetEndProportion;
		this.Lambda = new double[3];
		this.dataSplit = new ArrayList<Double>[3];
		this.weights = new double[3];
		this.firstNonZero = 0;

	}

	//endregion

	//region Properties

	public final ArrayList<Double> getData()
	{
		return data;
	}

	public final double GetAlpha(int i)
	{
		return Alpha[i];
	}

	public final double GetLambda(int i)
	{
		return Lambda[i];
	}

	public final MethodType getMethod()
	{
		return method;
	}
	public final void setMethod(MethodType value)
	{
		method = value;
	}

	public final int getForgetInitPeriods()
	{
		return forgetInitPeriods;
	}
	public final void setForgetInitPeriods(int value)
	{
		forgetInitPeriods = value;
	}

	public final int getForgetEndPeriods()
	{
		return forgetEndPeriods;
	}
	public final void setForgetEndPeriods(int value)
	{
		forgetEndPeriods = value;
	}

	public final double getForgetEndProportion()
	{
		return forgetEndProportion;
	}
	public final void setForgetEndProportion(double value)
	{
		forgetEndProportion = value;
	}

	public final double getFirstNonZero()
	{
		return firstNonZero;
	}
	public final void setFirstNonZero(double value)
	{
		firstNonZero = value;
	}

	//endregion

	//region Public Methods

	/**  Load data for calculation 
	 @param data data 
	*/
	public final void LoadData(ArrayList<Double> data)
	{
		this.data = data;
		boolean useMemory = (forgetInitPeriods > 0 && getForgetEndPeriods() > 0 && forgetEndPeriods < data.size() && forgetEndProportion > 0 && forgetEndProportion < 1);
		CalculateParameters(data, useMemory);
		CalculateMax();
	}

	/**  Probability of a number of events (or quantity) to happen in span period 
	 @param x number of events or quantity 
	 @param acum if probability is accumulated 
	 @return  the probability 
	*/
	public final double Probability(double x, boolean acum)
	{
		double prob = 0;
		double p = 0;
		//for (int i = 0; i < 3; i++) {
			switch (method)
			{
				//case MethodType.Poisson: p = dd.Poisson((int)x, Lambda[i], acum); break;
				case Hurdle:
					p = dd.HurdlePoisson((int)x, alpha, lambda, acum);
					break;
			}
			//prob += weights[i] * p;
		//}
		prob = p;
		if (Double.isInfinite(p))
		{
			prob = 0;
		}
		return prob;
	}

	private double Probability(double x)
	{
		if (x >= this.max)
		{
			return 1.00;
		}
		return Probability(x, true);
	}

	/**  Quantile for a certain probability (inverse function) 
	 @param p the probability of a quantity to occur 
	 @return  the quantile 
	*/
	public final double Quantile(double p)
	{
		if (p >= 0.99)
		{
			return this.max;
		}
		int it = 0;
		tangible.RefObject<Integer> tempRef_it = new tangible.RefObject<Integer>(it);
		double q1 = rs.MonotoneBisection(Probability, true, 0, max, p, epsilon, tempRef_it, maxIt);
	it = tempRef_it.argValue;
		double q2 = q1 + 1;
		if (Math.abs(Probability(q1) - p) <= Math.abs(Probability(q2) - p))
		{
			return Math.round(q1);
		}
		else
		{
			return Math.round(q2);
		}
	}

	/**  Get the histogram of the distribution 
	 @return  the histogram 
	*/
	public final Histogram GetHistogram()
	{
		Histogram hist = new Histogram(100);
		ArrayList<Double> freqs = new ArrayList<Double>();
		double p;
		double pTot = 0;
		for (int x = 0; x <= this.max; x++)
		{
			p = Probability(x, false);
			if (p <= 0.0005)
			{
				p = 0;
			}
			pTot += p;
			if (p > 0 && x > 0.5 && firstNonZero == 0)
			{
				firstNonZero = Math.round(pTot * Math.pow(10, 3)) / Math.pow(10, 3);
			}
			if (x == max && pTot < 1.00)
			{
				p += 1.0 - pTot;
			}
			freqs.add(p);
		}
		hist.LoadDist(freqs);
		return hist;
	}

	//endregion

	//region Private Methods

	private void CalculateN0t(ArrayList<Double> data, tangible.RefObject<Integer> n0, tangible.RefObject<Integer> t)
	{
		n0.argValue = 0;
		t.argValue = 0;
		for (double d : data)
		{
			if (d <= 0)
			{
				n0.argValue++;
			}
			else
			{
				t.argValue = t.argValue + (int)d;
			}
		}
	}

	private void CalculateParameters(ArrayList<Double> data, boolean useMemory)
	{
		if (useMemory)
		{
			CalculateParametersUseMemory(data, forgetInitPeriods, forgetEndPeriods, forgetEndProportion);
		}
		else
		{
			CalculateParameters(data);
		}
	}

	private void CalculateParameters(ArrayList<Double> data)
	{
		int n = data.size();
		int n0 = 0;
		int t = 0;
		tangible.RefObject<Integer> tempRef_n0 = new tangible.RefObject<Integer>(n0);
		tangible.RefObject<Integer> tempRef_t = new tangible.RefObject<Integer>(t);
		CalculateN0t(data, tempRef_n0, tempRef_t);
	t = tempRef_t.argValue;
	n0 = tempRef_n0.argValue;

		int maxLambda = 300;
		alpha = (double)n0 / (double)n;
		lambda = -1;
		tangible.RefObject<Double> tempRef_lambda = new tangible.RefObject<Double>(lambda);
		dd.MinimizeLikelihood(alpha, tempRef_lambda, maxLambda, n, n0, t);
	lambda = tempRef_lambda.argValue;
	}

	private void CalculateParametersUseMemory(ArrayList<Double> data, int forgetInitPeriods, int forgetEndPeriods, double forgetEndProportion)
	{
		double[] weights = func.CalcTempWeights(data.size(), forgetInitPeriods, forgetEndPeriods, forgetEndProportion);

		double n = 0;
		double n0 = 0;
		double t = 0;
		for (int i = 0;i < data.size();i++)
		{
			n += weights[i];
			if (data.get(i).compareTo(0) <= 0)
			{
				n0 += weights[i];
			}
			else
			{
				t = t + data.get(i) * weights[i];
			}
		}

		int maxLambda = 300;
		alpha = n0 / n;
		lambda = -1;
		int nInt = (int)Math.round(n);
		int n0Int = (int)Math.round(n0);
		int tInt = (int)Math.round(t);
		tangible.RefObject<Double> tempRef_lambda = new tangible.RefObject<Double>(lambda);
		dd.MinimizeLikelihood(alpha, tempRef_lambda, maxLambda, nInt, n0Int, tInt);
	lambda = tempRef_lambda.argValue;
	}

	private double CalculateMax()
	{
		this.max = func.Max(data);
		while (Probability(max) < 0.99)
		{
			if (max < 20)
			{
				max = max + 1;
			}
			else
			{
				max = max * 1.1;
			}
		}
		return max;
	}

	//endregion

	//region Obsolete

	private void SplitData(List<Double> rawData, int forgetInitPeriods, int forgetEndPeriods, double forgetEndProportion)
	{
		if (forgetInitPeriods == 0 && forgetEndPeriods == 0)
		{
			forgetInitPeriods = (int)(rawData.size() * 2.0 / 3.0);
			forgetEndPeriods = (int)(rawData.size() / 3.0);
			forgetEndProportion = 0.2;
		}
		dataSplit[0] = new ArrayList<Double>();
		dataSplit[1] = new ArrayList<Double>();
		dataSplit[2] = new ArrayList<Double>();
		for (int i = rawData.size() - 1; i >= 0; i--)
		{
			if (i > rawData.size() - forgetInitPeriods)
			{
				dataSplit[0].add(rawData.get(i));
			}
			else if (i > rawData.size() - forgetEndPeriods)
			{
				dataSplit[1].add(rawData.get(i));
			}
			else
			{
				dataSplit[2].add(rawData.get(i));
			}
		}

		weights[0] = 1 * dataSplit[0].size();
		weights[1] = ((1 - forgetEndProportion) / 2.0) * dataSplit[1].size();
		weights[2] = forgetEndProportion * dataSplit[2].size();
		weights = Normalize(weights);

	}

	/*
	private void CalculateLambdas(int span) {
	    alpha = new double[3];
	    switch (method) {
	        case MethodType.Poisson:
	            for (int i = 0; i < 3; i++) {
	                lambda[i] = CalculateLambda(dataSplit[i], span);
	            }
	            break;
	        case MethodType.Hurdle:
	            for (int i = 0; i < 3; i++) {
	                double alph = 0;
	                double lambd = 0;
	                CalculateAlphaLambda(dataSplit[i], span, ref alph, ref lambd);
	                alpha[i] = alph;
	                lambda[i] = lambd;
	            }
	            break;
	    }

	}
	*/

	private ArrayList<Double> Group(List<Double> rawData, int group)
	{
		max = 0;
		data = new ArrayList<Double>();
		int day = 0;
		double totPeriod = 0.0;
		double val;
		int nonZero = 0;
		for (int i = rawData.size() - 1; i >= 0; i--)
		{
			val = rawData.get(i);
			totPeriod += val;
			day++;
			if (day == group)
			{
				data.add(0, totPeriod);
				day = 0;
				if (totPeriod > 0)
				{
					nonZero++;
				}
				totPeriod = 0.0;
			}
		}
		if (day > 0)
		{
			data.add(0, totPeriod);
			if (totPeriod > 0)
			{
				nonZero++;
			}
			if (totPeriod > max)
			{
				max = totPeriod;
			}
		}
		loadFactor = (double)nonZero / (double)data.size();
		return data;
	}

	private double CalculateLambda(List<Double> data)
	{
		double lambda = stat.Mean(data);
		return lambda;
	}

	private double[] Normalize(List<Double> W)
	{
		double[] P = new double[W.size()];
		double sumW = 0;
		for (double w : W)
		{
			sumW += w;
		}
		for (int i = 0; i < W.size(); i++)
		{
			P[i] = W.get(i) / sumW;
		}
		return P;
	}

	private double CalculateLambda(List<Double> rawData, int span)
	{
		double sumLambda = 0;
		double lambda, maX;
		ArrayList<Double> rawSampData;
		ArrayList<Double> grSamp;
		for (int i = 0; i < maxIt; i++)
		{
			rawSampData = Resample(rawData);
			grSamp = Group(rawSampData, span);
			lambda = stat.Mean(grSamp);
			maX = func.Max(grSamp);
			if (maX > this.max)
			{
				this.max = maX;
			}
			sumLambda += lambda;
		}
		return sumLambda / (double)maxIt;
	}

	/*
	private void CalculateAlphaLambda(IList<double> rawData, int span, ref double alpha, ref double lambda) {
	    this.max = func.Max(rawData);

	    int n0 = 0;
	    int t = 0;
	    foreach (double d in grData) {
	        if (d == 0) { n0++; }
	        if (d > 0) { t = t + (int)d; }
	    }
	    int maxLambda = 300;
	    int n = rawData.Count;
	    dd.MinimizeLikelihood(ref alpha, ref lambda, maxLambda, n, n0, t);
	}
	*/
	/*
	private double CalculateMax() {
	    while (Probability(max) < 0.99) {
	        if (max < 20) { max = max + 1; } else { max = max * 1.1; }
	    }
	    return max;
	}
	*/
	private ArrayList<Double> Resample(List<Double> rawData, List<Double> W)
	{
		ArrayList<Double> resamp = new ArrayList<Double>();
		int n = rawData.size();
		int index;
		double u;

		while (resamp.size() < n)
		{
			index = rg.NextInt(0, n - 1);
			u = rg.NextDouble(0, 1);
			if (u <= W.get(index))
			{
				resamp.add(rawData.get(index));
			}
		}
		return resamp;
	}

	private ArrayList<Double> Resample(List<Double> rawData)
	{
		ArrayList<Double> resamp = new ArrayList<Double>();
		int index;
		while (resamp.size() < rawData.size())
		{
			index = rg.NextInt(0, rawData.size() - 1);
			resamp.add(rawData.get(index));
		}
		return resamp;
	}

	//endregion

	//region Enums

	public enum MethodType
	{
		Poisson,
		Hurdle;

		public int getValue()
		{
			return this.ordinal();
		}

		public static MethodType forValue(int value)
		{
			return values()[value];
		}
	}

	//endregion
}