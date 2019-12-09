package NumCalc;

import ClassicalStat.*;
import java.util.*;

//region Imports


//endregion


/**  Class for normal aproximation of convolutions 
*/
public class NormalConvolution extends ConvolutionCalc
{

	//region Fields

	private ArrayList<Double> data;
	private double n;
	private double convMean;
	private double convStDev;
	private Statistics stat;
	private NormalDist nd;

	//endregion

	//region Constructor

	/**  Constructor 
	*/
	public NormalConvolution()
	{
		data = new ArrayList<Double>();
		stat = new Statistics();
		nd = new NormalDist();
	}

	//endregion

	//region Properties

	/**  Convolution mean 
	*/
	public final double getConvMean()
	{
		return convMean;
	}

	/**  Convolution standard deviation 
	*/
	public final double getConvStDev()
	{
		return convStDev;
	}

	//endregion


	//region ConvolutionCalc interface implementation

	/**  Load data for calculation 
	 @param data data from time series 
	 @param n number (real) of convolutions 
	*/
	public final void LoadData(ArrayList<Double> data, double n)
	{
		this.data = data;
		this.n = n;
		this.convMean = stat.Mean(data) * n;
		this.convStDev = stat.StDev(data) * Math.sqrt(n);
	}

	/**  Load histogram for calculation 
	 @param hist histogram for calculation 
	 @param n number of convolutions 
	*/
	public final void LoadHistogram(Histogram hist, double n)
	{
		this.n = n;
		this.convMean = hist.Mean * n;
		this.convStDev = hist.StDev * Math.sqrt(n);
	}

	/**  Acumulated probability for a certain value 
	 @param x the value 
	 @return  the probability 
	*/
	public final double ProbabilityAcum(double x)
	{
		if (convStDev == 0)
		{
			if (x < convMean)
			{
				return 0;
			}
			else
			{
				return 1;
			}
		}
		return nd.pNorm(x, convMean, convStDev);
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
			throw new RuntimeException(Properties.Strings.probability_must_be_between_0_and_1);
		}
		return nd.qNorm(p, convMean, convStDev);
	}

	/**  if calculation is valid 
	 @return  true if it is valid, false if not 
	*/
	public final boolean Valid()
	{
		return true;
	}

	//endregion
}