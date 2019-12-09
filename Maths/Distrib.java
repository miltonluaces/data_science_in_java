package NumCalc;

//region Imports

import ClassicalStat.*;
import java.util.*;

//endregion


/**  Class for Distribution approximation and statistic calculations based in  and trimmed/winsored normal approximation and/or kernels 
*/
public class Distrib {

	//region Fields

	private int maxForKernels;
	private double s;
	private int maxIterations;
	private int maxClasses;
	private Statistics stat;
	private NormalDist normalDist;
	private DensityEstim dEst;
	private boolean byKernels;
	private boolean windsor;
	private double cutoff;
	private double mean;
	private double sd;

	//endregion

	//region Constructor

	/**  Constructor 
	 @param maxForKernels maximum of cases for kernel calculations 
	 @param s standard deviation for kernels 
	 @param maxIterations maximum of iterations for kernels 
	 @param maxClasses maximum of classes for kernels 
	 @param windsor if windorized sample (else trimmed sample) for larger data sets 
	 @param cutoff cut-off for trim and winsdor samples for larger data sets 
	*/
	public Distrib(int maxForKernels, double s, int maxIterations, int maxClasses, boolean windsor, double cutoff)
	{
		this.maxForKernels = maxForKernels;
		this.s = s;
		this.maxIterations = maxIterations;
		this.maxClasses = maxClasses;
		this.byKernels = false;
		this.windsor = false;
		this.cutoff = cutoff;
		this.mean = -1;
		this.sd = -1;
		this.stat = new Statistics();
		this.normalDist = new NormalDist();
		this.dEst = new DensityEstim(s, maxIterations, maxClasses);
	}

	//endregion

	//region Public Methods

	//region Properties

	/**  maximum number of cases for using kernels 
	*/
	public final int getMaxForKernels()
	{
		return maxForKernels;
	}
	public final void setMaxForKernels(int value)
	{
		maxForKernels = value;
	}

	/**  standard deviation for kernels 
	*/
	public final double getS()
	{
		return s;
	}
	public final void setS(double value)
	{
		s = value;
	}

	/**  maximum number of iterations for kernels 
	*/
	public final int getMaxIterations()
	{
		return maxIterations;
	}
	public final void setMaxIterations(int value)
	{
		maxIterations = value;
	}

	/**  maximum number of classes for kernels 
	*/
	public final int getMaxClasses()
	{
		return maxClasses;
	}
	public final void setMaxClasses(int value)
	{
		maxClasses = value;
	}

	/**  if the sample should be winsored (if not, it should be trimmed) 
	*/
	public final boolean getWindsor()
	{
		return windsor;
	}
	public final void setWindsor(boolean value)
	{
		windsor = value;
	}

	/**  cut-off for both the trimmed or winsored statistics 
	*/
	public final double getCutoff()
	{
		return cutoff;
	}
	public final void setCutoff(double value)
	{
		cutoff = value;
	}

	/**  if the calculation was done by kernels (read-only) 
	*/
	public final boolean getByKernels()
	{
		return byKernels;
	}

	/**  robust mean calculated (read-only) 
	*/
	public final double getMean()
	{
		return mean;
	}

	/**  robust standard deviation calculated (read-only) 
	*/
	public final double getSd()
	{
		return sd;
	}

	//endregion

	//region Load data

	/**  Load data for ditribution calculations 
	 @param data data for load 
	 @param sort if data should be sorted or not 
	*/
	public final void LoadData(ArrayList<Double> data, boolean sort)
	{
		byKernels = (data.size() < maxForKernels);
		if (byKernels)
		{
			LoadKernels(data);
			dEst.SetPercentiles(false);
		}
		else
		{
			ArrayList<Double> robData = new ArrayList<Double>();
			if (windsor)
			{
				robData = stat.Windsor(data, cutoff, sort);
			}
			else
			{
				robData = stat.Trim(data, cutoff, sort);
			}
			mean = stat.Mean(robData);
			sd = stat.StDev(robData);
		}
	}

	//endregion 

	//region Calculate statistics 

	/**  Probability of a certain value 
	 @param val the value 
	 @return  the calculated probability 
	*/
	public final double Probability(double val)
	{
		if (byKernels)
		{
			return dEst.Probability(val);
		}
		else
		{
			return normalDist.pNorm(val, mean, sd);
		}
	}

	/**  Quantile of a certain probability 
	 @param p the probability 
	 @return  the calcualted quantile 
	*/
	public final double Quantile(double p)
	{
		if (p > 0.9999)
		{
			return mean + 4 * sd;
		}
		if (byKernels)
		{
			return dEst.CalculatePercentile(p * 100);
		}
		else
		{
			return normalDist.qNorm(p, mean, sd);
		}
	}

	//endregion

	//endregion

	//region Private Methods

	private void LoadKernels(ArrayList<Double> data)
	{
		int max = Integer.MIN_VALUE;
		for (int val : data)
		{
			if (val > max)
			{
				max = val;
			}
		}
		double[] meanFrec = new double[max + 1];
		for (int val : data)
		{
			meanFrec[val] = meanFrec[val] + 1;
		}
		dEst = new DensityEstim(s, maxIterations, maxClasses);
		dEst.LoadDist(new ArrayList<Double>(meanFrec));
	}

	//endregion
}