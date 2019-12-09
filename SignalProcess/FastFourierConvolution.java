package SignalProcess;

import ClassicalStat.*;
import NumCalc.*;
import Fourier.Properties.*;
import java.util.*;

//region Imports


//endregion


public class FastFourierConvolution extends ConvolutionCalc
{

	//region Fields

	private ArrayList<Double> probs;
	private ArrayList<Double> acumProbs;
	private FFT fft;
	private Histogram histogram;
	private double n;
	private int maxClasses;
	private int lag;

	//endregion

	//region Constructor

	/**  Constructor 
	 @param maxClasses maximum number of classes 
	*/
	public FastFourierConvolution(int maxClasses)
	{
		this.fft = new FFT();
		this.maxClasses = maxClasses;
	}

	//endregion

	//region Public Methods

	/**  Load data for calculation 
	 @param data data from time series 
	 @param n number of convolutions 
	*/
	public final void LoadData(ArrayList<Double> data, double n)
	{
		histogram.Clear();
		histogram.LoadData(data);
		this.n = n;
	}

	/**  Load histogram for calculation 
	 @param hist histogram for calculation 
	 @param n number of convolutions 
	*/
	public final void LoadHistogram(Histogram hist, double n)
	{
		histogram.Clear();
		this.histogram = hist;
		this.n = n;
	}

	/**  Acumulated probability for a certain value 
	 @param x the value 
	 @return  the probability 
	*/
	public final double ProbabilityAcum(double x)
	{
		if (probs == null)
		{
			LinearConvolve(histogram, (int)n);
		}
		int key = (int)Math.ceil(((x / histogram.Width - histogram.Min + histogram.Mean)));
		if (key < 0 || key > probs.size() - 1)
		{
			return 0;
		}
		return acumProbs.get(key);
	}

	/**  Quantile for a certain probability 
	 @param p the probability 
	 @return  the value 
	*/
	public final double Quantile(double p)
	{
		if (p > 1)
		{
			throw new RuntimeException(Strings.Error_Probabilities_must_be_between_0_and_1);
		}
		if (probs == null)
		{
			LinearConvolve(histogram, (int)n);
		}

		for (int i = 0; i < acumProbs.size(); i++)
		{
			if (acumProbs.get(i).compareTo(p) >= 0)
			{
				if (i > 0 && Math.abs(acumProbs.get(i) - p) > Math.abs(acumProbs.get(i + 1) - p))
				{
					return (int)Math.round(i + 1 - histogram.Mean);
				}
				else
				{
					return (int)Math.round(i - histogram.Mean);
				}
			}
		}
		return -1;
	}

	/**  if calculation is valid 
	 @return  true if it is valid, false if not 
	*/
	public final boolean Valid()
	{
		return true;
	}

	//endregion

	//region Private Methods

	private double ProbabilityFft(int n, double value)
	{
		if (probs == null)
		{
			LinearConvolve(histogram, n);
		}
		//int key = (int)(value / histogram.Width- histogram.Min) + lag;
		int key = (int)Math.ceil(((value / histogram.Width - histogram.Min) + histogram.Mean));
		if (key < 0 || key > probs.size() - 1)
		{
			return 0;
		}
		return probs.get(key);
	}


	private void LinearConvolve(Histogram hist, int n)
	{
		ArrayList<Double> completeValues = new ArrayList<Double>();
		for (int i = 0; i < maxClasses; i++)
		{
			if (hist.GetValue(i) > hist.Max)
			{
				break;
			}
			if (hist.Freqs.ContainsKey(i))
			{
				completeValues.add(hist.Freqs[i]);
			}
			else
			{
				completeValues.add(0.0);
			}
		}
		PadRight(completeValues);

		Complex[] comp = new Complex[completeValues.size()];
		for (int i = 0; i < completeValues.size(); i++)
		{
			comp[i] = new Complex(completeValues.get(i), 0);
		}
		Complex[] complexRes = fft.LinearConvolve(comp, n);
		ArrayList<Double> realRes = new ArrayList<Double>();
		double tot = 0;
		double val;
		for (int i = 0; i < complexRes.length; i++)
		{
			val = complexRes[i].Real;
			if (val > 0)
			{
				realRes.add(val);
				tot += val;
			}
			else
			{
				realRes.add(0);
			}

		}
		probs = new ArrayList<Double>();
		acumProbs = new ArrayList<Double>();
		double prob;
		double acum = 0;
		for (int i = 0; i < complexRes.length; i++)
		{
			prob = realRes.get(i) / tot;
			probs.add(realRes.get(i) / tot);
			acum += prob;
			acumProbs.add(acum);
		}
		lag = (int)Math.round(hist.Mean);
	}

	private void PadRight(ArrayList<Double> values)
	{
		int count = values.size();
		int quadPow = 2;
		while (count > quadPow)
		{
			quadPow *= 2;
		}
		if (quadPow == count)
		{
			return;
		}
		int addCount = quadPow - count;
		for (int i = 0; i < addCount; i++)
		{
			values.add(0);
		}
	}

	//endregion
}