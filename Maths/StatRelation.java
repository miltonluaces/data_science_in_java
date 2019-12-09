package NumCalc;

//region Imports

import ClassicalStat.*;
import java.util.*;

//endregion


public class StatRelation
{

	//region Fields

	private Statistics stat;
	private NormalDist nd;
	private DensityEstim de;
	private double s;
	private int maxIterations;
	private int maxClasses;
	private Splines sp;
	private TreeMap<Integer, ArrayList<Double>> relation;
	private int factor;
	private MethodType method = MethodType.values()[0];
	private boolean interpolate;
	private ANOVA anova;

	//endregion

	//region Constructor

	public StatRelation()
	{
		factor = 100;
		method = MethodType.normal;
		stat = new Statistics();
		nd = new NormalDist();
		de = new DensityEstim(s, maxIterations, maxClasses);
		sp = new Splines();
		relation = new TreeMap<Integer, ArrayList<Double>>();
		factor = 100;
		interpolate = false;
		anova = new ANOVA();

	}

	//endregion

	//region Properties

	public final int getFactor()
	{
		return factor;
	}
	public final void setFactor(int value)
	{
		factor = value;
	}
	public final MethodType getMethod()
	{
		return method;
	}
	public final void setMethod(MethodType value)
	{
		method = value;
	}

	public final boolean getInterpolate()
	{
		return interpolate;
	}
	public final void setInterpolate(boolean value)
	{
		interpolate = value;
	}

	public final double getS()
	{
		return s;
	}
	public final void setS(double value)
	{
		s = value;
	}

	public final int getMaxIterations()
	{
		return maxIterations;
	}
	public final void setMaxIterations(int value)
	{
		maxIterations = value;
	}

	public final int getMaxClasses()
	{
		return maxClasses;
	}
	public final void setMaxClasses(int value)
	{
		maxClasses = value;
	}

	//endregion

	//region Public Methods

	public final void LoadData(ArrayList<Double> X, ArrayList<Double> Y)
	{
		if (X.size() != Y.size())
		{
			int minCount = Math.min(X.size(), Y.size());
			while (X.size() > minCount)
			{
				X.remove(0);
			}
			while (Y.size() > minCount)
			{
				Y.remove(0);
			}
		}
		for (int i = 0;i < X.size();i++)
		{
			int key = (int)(X.get(i) * factor);
			if (!relation.containsKey(key))
			{
				relation.put(key, new ArrayList<Double>());
			}
			relation.get(key).add(Y.get(i));
		}
	}

	public final double Anova()
	{
		List<List<Double>> groups = new ArrayList<List<Double>>();
		for (ArrayList<Double> group : relation.values())
		{
			groups.add(group);
		}
		double pValue = anova.Calculate(groups);
		return pValue;
	}

	public final ArrayList<Double> GetX()
	{
		ArrayList<Double> X = new ArrayList<Double>();
		for (int key : relation.keySet())
		{
			X.add((double)key / factor);
		}
		return X;
	}

	public final ArrayList<Double> GetY(double x)
	{
		int key = GetXKey(x);
		if (!relation.containsKey(key))
		{
			throw new RuntimeException("Error. Key not found");
		}
		return relation.get(key);
	}

	public final double Estimate(double x, double p)
	{
		double q = -1;
		int xKey = GetXKey(x);
		int lowKey = -1;
		int uppKey = -1;
		for (int key : relation.keySet())
		{
			if (lowKey == -1)
			{
				lowKey = key;
			}
			uppKey = key;
			if (uppKey > xKey)
			{
				break;
			}
		}
		if (interpolate)
		{
			double lowQ = GetQ(relation.get(lowKey), p);
			double uppQ = GetQ(relation.get(uppKey), p);
			return sp.Interpolar(xKey, lowKey, lowQ, uppKey, uppQ);
		}
		else
		{
			if (xKey - lowKey < uppKey - xKey)
			{
				q = GetQ(relation.get(lowKey), p);
			}
			else
			{
				q = GetQ(relation.get(uppKey), p);
			}
		}
		return q;
	}

	//endregion

	//region Private Methods

	private double GetQ(ArrayList<Double> data, double p)
	{
		switch (method)
		{
			case normal:
				return GetNormalQ(data, p);
			case empirical:
				return GetKernelQ(data, p);
			default:
				throw new RuntimeException("Error. Method not found");
		}
	}

	private double GetKernelQ(ArrayList<Double> data, double p)
	{
		Histogram hist = new Histogram(maxClasses);
		hist.LoadData(data);
		de.LoadHist(hist);
		return de.CalculatePercentile(p);
	}

	private double GetNormalQ(ArrayList<Double> data, double p)
	{
		double mean = stat.Mean(data);
		double sd = stat.StDev(data);
		return nd.qNorm(p, mean, sd);
	}

	private int GetXKey(double value)
	{
		return (int)(value * factor);
	}
	private double GetXValue(int key)
	{
		return key / factor;
	}

	//endregion

	//region Public Enums

	public enum MethodType
	{
		normal,
		empirical;

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