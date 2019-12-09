package NumCalc;

import ClassicalStat.*;
import java.util.*;

//region Imports


//endregion


/**  Class for convolution calculation based on previous normalizing of data 
*/
public class NormalizedConvolution extends ConvolutionCalc
{

	//region Fields

	private NormalDist norm;
	private HypoTest hypo;
	private Statistics stat;
	private double mean;
	private double stDev;
	private double convMean;
	private double convStDev;
	private double lambda;
	private double pValue;
	private double minMaxLik;
	private double n;
	private double maxLik;
	private double lambdaInc;

	//endregion

	//region Constructor

	/**  Constructor 
	 @param lambdaInc increments to apply on lambda from -3 to 3 (Box-Cox family) 
	 @param minMaxLik minimum value for maximum likelihood 
	*/
	public NormalizedConvolution(double lambdaInc, double minMaxLik)
	{
		norm = new NormalDist();
		stat = new Statistics();
		hypo = new HypoTest(0.05, 10);
		this.lambdaInc = lambdaInc;
		pValue = -1;
		maxLik = -1;
		this.minMaxLik = minMaxLik;
	}

	//endregion

	//region Properties

	/**  Maximum likelihood value 
	*/
	public final double getMaxLik()
	{
		return maxLik;
	}

	/**  Selected lambda value 
	*/
	public final double getLambda()
	{
		return lambda;
	}

	//endregion

	//region ConvolutionCalc interface implementation

	/**  Load data for calculation 
	 @param data data from time series 
	 @param n number of convolutions 
	*/
	public final void LoadData(ArrayList<Double> data, double n)
	{
		this.n = n;
		maxLik = -1;
		tangible.RefObject<Double> tempRef_maxLik = new tangible.RefObject<Double>(maxLik);
		lambda = hypo.GetBoxCoxLambda(data, lambdaInc, tempRef_maxLik);
	maxLik = tempRef_maxLik.argValue;
		ArrayList<Double> norm = hypo.GetBoxCoxTransf(data, lambda);
		mean = stat.Mean(norm);
		stDev = stat.StDev(norm);
		convMean = n * mean;
		convStDev = Math.sqrt(n) * stDev;
	}

	/**  Load histogram for calculation 
	 @param hist histogram for calculation 
	 @param n number of convolutions 
	*/
	public final void LoadHistogram(Histogram hist, double n)
	{
		LoadData(hist.GetRawData(), n);
	}

	/**  Acumulated probability for a certain value 
	 @param x the value 
	 @return  the probability 
	*/
	public final double ProbabilityAcum(double x)
	{
		double xTr = hypo.BoxCoxTransf(x, lambda);
		if (convStDev == 0)
		{
			if (xTr < convMean)
			{
				return 0;
			}
			else
			{
				return 1;
			}
		}
		return norm.pNorm(xTr, convMean, convStDev);
	}

	/**  Quantile for a certain probability 
	 @param p the probability 
	 @return  the value 
	*/
	public final double Quantile(double p)
	{
		if (convStDev == 0)
		{
			return convMean;
		}
		if (p < 0 || p > 1)
		{
			throw new RuntimeException(Properties.Strings.Error_ProbabilityMustBeBetween0_1);
		}
		double valueTr = norm.qNorm(p, convMean, convStDev);
		return hypo.BoxCoxTransfInv(valueTr, lambda);

	}

	/**  if calculation is valid 
	 @return  true if it is valid, false if not 
	*/
	public final boolean Valid()
	{
		//return (maxLik >= minMaxLik);  
		return false;
	}

	//endregion

}