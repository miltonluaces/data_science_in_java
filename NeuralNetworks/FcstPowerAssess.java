package NeuralNetworks;

import NumCalc.*;
import ClassicalStat.*;
import java.util.*;

//region Imports


//endregion


/**  Class for power assess calculations 
*/
public class FcstPowerAssess
{

	//region Fields

	private ArrayList<Double> hist;
	private ArrayList<Double> fcst;

	private double trendTreshold;

	private int movingWindow;
	private int iterations;
	private double serviceLevel;

	private int count;
	private double hits;
	private double neut;
	private double miss;

	private TDNN tdnn;
	private Functions func;
	private Polynomic poly;
	private Statistics stat;
	private NormalDist normalDist;

	private ArrayList<Double> assessList;

	private double valid;
	private double power;

	private double meanDiffStock;
	private double meanSLevel;
	private double meanBackorder;

	private FcstMethodType fcstMethod = FcstMethodType.values()[0];
	private AssessMethodType assessMethod = AssessMethodType.values()[0];

	//endregion

	//region Constructor

	/**  Constructor 
	*/
	public FcstPowerAssess()
	{
		this.func = new Functions();
		this.poly = new Polynomic();
		this.stat = new Statistics();
		this.normalDist = new NormalDist();
		this.trendTreshold = 0.2;
		this.movingWindow = 7;
		this.iterations = 50;
		this.serviceLevel = 0.95;

		this.fcstMethod = FcstMethodType.LtFcNaive;
		this.assessMethod = AssessMethodType.StError;
	}

	//endregion

	//region Properties

	/**  Historic grouped series 
	*/
	public final ArrayList<Double> getHist()
	{
		return hist;
	}

	/**  Forecast grouped series 
	*/
	public final ArrayList<Double> getFcst()
	{
		return fcst;
	}

	/**  Threshold used for tend assessment 
	*/
	public final double getTrendThreshold()
	{
		return trendTreshold;
	}
	public final void setTrendThreshold(double value)
	{
		trendTreshold = value;
	}

	/**  Moving window for moving average and other algorithms 
	*/
	public final int getMovingWindow()
	{
		return movingWindow;
	}
	public final void setMovingWindow(int value)
	{
		movingWindow = value;
	}

	/**  Number of iterations in case of iterative algorithms 
	*/
	public final int getIterations()
	{
		return iterations;
	}
	public final void setIterations(int value)
	{
		iterations = value;
	}

	/**  Number of data 
	*/
	public final double getCount()
	{
		return count;
	}

	/**  Service level for comparison of methods 
	*/
	public final double getServiceLevel()
	{
		return serviceLevel;
	}
	public final void setServiceLevel(double value)
	{
		serviceLevel = value;
	}

	/**  Hit cases for trending assessment 
	*/
	public final double getHits()
	{
		return hits;
	}

	/**  Neutral cases for trending assessment 
	*/
	public final double getNeut()
	{
		return neut;
	}

	/**  Miss cases for trending assessment 
	*/
	public final double getMiss()
	{
		return miss;
	}

	/**  Percentage of valid periods for assessment (hits or miss, not  neutral) 
	*/
	public final double getValid()
	{
		return valid;
	}

	/**  Percentage of forecasting power 
	*/
	public final double getPower()
	{
		return power;
	}

	/**  Mean stock 
	*/
	public final double getMeanDiffStock()
	{
		return meanDiffStock;
	}

	/**  Mean service level 
	*/
	public final double getMeanSLevel()
	{
		return meanSLevel;
	}

	/**  Mean backorder 
	*/
	public final double getMeanBackorder()
	{
		return meanBackorder;
	}

	/**  Forecasting method for calculation 
	*/
	public final FcstMethodType getFcstMethod()
	{
		return fcstMethod;
	}
	public final void setFcstMethod(FcstMethodType value)
	{
		fcstMethod = value;
	}

	/**  Assessment method for calculation 
	*/
	public final AssessMethodType getAssessMethod()
	{
		return assessMethod;
	}
	public final void setAssessMethod(AssessMethodType value)
	{
		assessMethod = value;
	}

	/**  List of assessment values 
	*/
	public final ArrayList<Double> getAssessList()
	{
		return assessList;
	}

	//endregion

	//region Public Methods

	//region Load Data

	/**  Load data for calculation 
	 @param hist historic time series (grouped) 
	 @param fcst forcast time series (grouped) 
	*/
	public final void LoadData(ArrayList<Double> hist, ArrayList<Double> fcst)
	{
		if (hist.size() != fcst.size())
		{
			throw new RuntimeException("Error");
		}
		this.count = hist.size();
		this.hist = hist;
		this.fcst = fcst;
	}

	/**  Load data for calculation 
	 @param hist historic time series (grouped) 
	*/
	public final void LoadData(ArrayList<Double> hist)
	{
		this.hist = hist;
		this.count = hist.size();
		switch (fcstMethod)
		{
			case Naive:
				this.fcst = CalcNaiveFcst(hist);
				break;
			case MovAvg:
				this.fcst = CalcMAvgFcst(hist, movingWindow, false);
				break;
			case Regress:
				this.fcst = CalcRegressFcst(hist, movingWindow);
				break;
			case TDNN:
				this.fcst = CalcTDNNFcst(hist, movingWindow, iterations);
				break;
			case LtFcNaive:
				this.fcst = CalcLtNaiveFcst(hist, serviceLevel);
				break;
			case LtFcMovAvg:
				this.fcst = CalcLtMAvgFcst(hist, movingWindow, false, serviceLevel);
				break;
		}
	}

	//endregion

	//region Calculate

	/**  Assess using assessment method 
	*/
	public final void Assess()
	{
		switch (assessMethod)
		{
			case Trend:
				TrendAssess();
				break;
			case StError:
				StErrorAssess();
				break;
			case Metrics:
				MetricsAssess();
				break;
		}

	}

	//endregion

	//endregion

	//region Private Methods

	//region Assessment methods

	private void TrendAssess()
	{
		hits = 0;
		neut = 0;
		miss = 0;
		assessList = new ArrayList<Double>();

		ArrayList<Double> trendHist = GetTrend(hist, true);
		ArrayList<Double> trendFcst = GetTrend(fcst, false);

		for (int i = 0; i < hist.size(); i++)
		{
			if (trendFcst.get(i).equals(-1) || (trendHist.get(i).equals(0) && trendFcst.get(i).equals(0)))
			{
				assessList.add(0.0);
				neut++;
				continue;
			}
			if ((int)Math.signum(trendHist.get(i)) == (int)Math.signum(trendFcst.get(i)))
			{
				assessList.add(1.0);
				hits++;
			}
			else
			{
				assessList.add(-1.0);
				miss++;
			}
		}
		valid = (double)(hits + miss) / (double)count;
		power = (double)(hits) / (hits + miss);
	}

	private void StErrorAssess()
	{
		assessList = new ArrayList<Double>();

		double sumHist = 0;
		double sumErr = 0;
		double err;
		for (int i = 0; i < hist.size() - 1; i++)
		{
			if (fcst.get(i).equals(-1))
			{
				continue;
			}
			err = Math.abs(hist.get(i + 1) - fcst.get(i));
			assessList.add(err);
			sumHist += hist.get(i + 1);
			sumErr += err;
		}
		power = 1.0 - sumErr / sumHist;
	}

	private void MetricsAssess()
	{
		double sumSL = 0;
		double sumBack = 0;
		double sumDiffStock = 0;
		for (int i = 0; i < hist.size() - 1; i++)
		{
			if (fcst.get(i).equals(-1))
			{
				continue;
			}
			if (fcst.get(i).compareTo(hist.get(i + 1)) >= 0)
			{
				sumSL++;
				sumDiffStock += fcst.get(i) - hist.get(i + 1);
			}
			else
			{
				sumBack += hist.get(i + 1) - fcst.get(i);
			}

		}
		double n = (double)(hist.size() - 1);
		meanDiffStock = sumDiffStock / n;
		meanSLevel = sumSL / n;
		meanBackorder = sumBack / n;
	}

	//endregion 

	//region Basic Forecast calculations

	private ArrayList<Double> CalcNaiveFcst(ArrayList<Double> hist)
	{
		return hist;
	}

	private ArrayList<Double> CalcMAvgFcst(ArrayList<Double> hist, int mw, boolean centered)
	{
		ArrayList<Double> fcst = func.MovingAverage(hist, mw, centered);
		return fcst;
	}

	private ArrayList<Double> CalcLtNaiveFcst(ArrayList<Double> hist, double serviceLevel)
	{
		double mean = stat.Mean(hist);
		double nStDevs = normalDist.qNorm(serviceLevel / 100.0, 0, 1);
		ArrayList<Double> ltFcst = new ArrayList<Double>();
		double stDev = stat.StDev(hist);
		for (int i = 0; i < hist.size(); i++)
		{
			ltFcst.add(mean + nStDevs * stDev);
		}
		return ltFcst;
	}

	private ArrayList<Double> CalcLtMAvgFcst(ArrayList<Double> hist, int mw, boolean centered, double serviceLevel)
	{
		double nStDevs = normalDist.qNorm(serviceLevel, 0, 1);
		double stDev = stat.StDev(hist);
		ArrayList<Double> fcst = func.MovingAverage(hist, mw, centered);
		ArrayList<Double> ltFcst = new ArrayList<Double>();
		for (int i = 0; i < fcst.size(); i++)
		{
			ltFcst.add(fcst.get(i) + nStDevs * stDev);
		}
		return ltFcst;
	}

	private ArrayList<Double> CalcRegressFcst(ArrayList<Double> hist, int mw)
	{
		ArrayList<Double> X = new ArrayList<Double>();
		for (int i = 0; i < hist.size(); i++)
		{
			X.add((double)i);
		}
		Result res = poly.LSRegression(hist);
		ArrayList<Double> coeffs = new ArrayList<Double>();
		coeffs.add(res.c0);
		coeffs.add(res.c1);
		ArrayList<Double> fcst = poly.GetRegValues(X, coeffs);
		return fcst;
	}

	private ArrayList<Double> CalcTDNNFcst(ArrayList<Double> hist, int mw, int epochs)
	{
		TDNN tdnn = new TDNN(mw);
		tdnn.LoadSerie(hist, 0);
		ArrayList<Double> fcst = tdnn.ModelFromHistory(epochs);
		return fcst;
	}

	//endregion

	//region Inverse moving window calculations

	private ArrayList<Double> CalcFcstPower(ArrayList<Double> hist, FcstMethodType fcstMethod, int minMw, int maxMw, boolean centered)
	{
		if (fcstMethod != FcstMethodType.MovAvg || fcstMethod != FcstMethodType.TDNN)
		{
			throw new RuntimeException("Error. Not implemented for this method");
		}

		ArrayList<Double> fcstPower = new ArrayList<Double>();
		this.hist = hist;
		this.count = hist.size();
		for (int mw = maxMw; mw >= minMw; mw--)
		{
			this.movingWindow = mw;
			LoadData(hist);
			TrendAssess();
			fcstPower.add(this.getPower());
		}
		return fcstPower;
	}

	//endregion

	//region Get Trend

	private ArrayList<Double> GetTrend(ArrayList<Double> serie, boolean next)
	{
		ArrayList<Double> trendSerie = new ArrayList<Double>();
		trendSerie.add(0.0);
		double tr;
		for (int i = 1; i < serie.size() - 1; i++)
		{
			if (next)
			{
				if (serie.get(i).equals(0))
				{
					if (serie.get(i + 1).compareTo(0.0) > 0)
					{
						tr = 1;
					}
					else if (serie.get(i + 1).compareTo(0.0) < 0)
					{
						tr = -1;
					}
					else
					{
						tr = 0;
					}
				}
				else
				{
					tr = (serie.get(i + 1) - serie.get(i)) / serie.get(i);
				}
			}
			else
			{
				if (serie.get(i - 1).equals(0))
				{
					if (serie.get(i).compareTo(0.0) > 0)
					{
						tr = 1;
					}
					else if (serie.get(i).compareTo(0.0) < 0)
					{
						tr = -1;
					}
					else
					{
						tr = 0;
					}
				}
				else
				{
					tr = (serie.get(i) - serie.get(i - 1)) / serie.get(i - 1);
				}
			}

			if (Math.abs(tr) < trendTreshold)
			{
				tr = 0;
			}
			trendSerie.add(tr);
		}
		trendSerie.add(0.0);
		return trendSerie;
	}

	//endregion

	//endregion

	//region Enums

	public enum FcstMethodType
	{
		Naive,
		MovAvg,
		Regress,
		TDNN,
		LtFcNaive,
		LtFcMovAvg;

		public int getValue()
		{
			return this.ordinal();
		}

		public static FcstMethodType forValue(int value)
		{
			return values()[value];
		}
	}
	public enum AssessMethodType
	{
		Trend,
		StError,
		Metrics;

		public int getValue()
		{
			return this.ordinal();
		}

		public static AssessMethodType forValue(int value)
		{
			return values()[value];
		}
	}

	//endregion
}