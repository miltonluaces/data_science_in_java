package NumCalc;

import ClassicalStat.*;
import java.util.*;

//region Imports


//endregion


/**  Class for trend determination (obsolete) 
*/
public class Trend
{


	//region Fields

	private ArrayList<Double> resids;
	private int minLeadTimeObs;
	private double trendUpThreshold;
	private double trendDownThreshold;
	private double varCoeffThreshold = 1;
	private double outlierThreshold;
	private Functions func;
	private Polynomic poly;
	private Statistics stat;

	//endregion

	//region Constructor

	/**  Constructor 
	 @param minLeadTimeObs minimum of observations 
	 @param trendUpThreshold threshold for ascending trend 
	 @param trendDownThreshold threshold for descending trend 
	 @param outlierThreshold threshold for outlier filtering 
	*/
	public Trend(int minLeadTimeObs, double trendUpThreshold, double trendDownThreshold, double outlierThreshold)
	{
		this.minLeadTimeObs = minLeadTimeObs;
		this.trendUpThreshold = trendUpThreshold;
		this.trendDownThreshold = trendDownThreshold;
		this.outlierThreshold = outlierThreshold;
		func = new Functions();
		poly = new Polynomic();
		stat = new Statistics();
	}

	//endregion

	//region Properties

	/**  Bias 
	*/
	public final ArrayList<Double> getResids()
	{
		return resids;
	}
	public final void setResids(ArrayList<Double> value)
	{
		resids = value;
	}

	/**  Minimum of observations 
	*/
	public final int getMinLeadTimeObs()
	{
		return minLeadTimeObs;
	}
	public final void setMinLeadTimeObs(int value)
	{
		minLeadTimeObs = value;
	}

	/**   Threshold for ascending trend 
	*/
	public final double getTrendUpThreshold()
	{
		return trendUpThreshold;
	}
	public final void setTrendUpThreshold(double value)
	{
		trendUpThreshold = value;
	}

	/**   Threshold for descending trend 
	*/
	public final double getTrendDownThreshold()
	{
		return trendDownThreshold;
	}
	public final void setTrendDownThreshold(double value)
	{
		trendDownThreshold = value;
	}

	//endregion

	//region Methods

	//region Function Trend

	/**  Trunk a time series 
	 @param serie the time series 
	 @param firstIndex first index for trunk 
	 @return  trunked serie 
	*/
	public final ArrayList<Double> Trunk(ArrayList<Double> serie, int firstIndex)
	{
		ArrayList<Double> trunk = new ArrayList<Double>();
		for (int i = firstIndex;i < serie.size();i++)
		{
			trunk.add(serie.get(i));
		}
		return trunk;
	}

	/**  Move first period due to descending trend 
	 @param serie the time series 
	 @return  the time series with moved first period 
	*/
	public final ArrayList<Double> MoveIniByDownTrend(ArrayList<Double> serie)
	{
		if (serie.size() <= minLeadTimeObs)
		{
			return serie;
		}
		ArrayList<Double> drv2 = poly.DeriveSerie(serie, 2);
		ArrayList<Double> trunk = serie;
		for (int index : poly.ZerosNegPosIndexes(drv2))
		{
			if (serie.size() - index < minLeadTimeObs)
			{
				return trunk;
			}
			trunk = new ArrayList<Double>();
			for (int i = index;i < serie.size();i++)
			{
				trunk.add(serie.get(i));
			}
			if (GetFunctionTrend(trunk) == 0)
			{
				return trunk;
			}
		}
		//if there are no more candidates
		int minIndex = 0;
		int maxIndex = trunk.size() - minLeadTimeObs;
		ArrayList<Double> tr = Trunk(trunk, maxIndex);

		//if it can be stationarized
		minIndex = GetMinIndex(trunk, 0, trunk.size() - 1, false, 0);
		trunk = Trunk(trunk, minIndex);
		return trunk;
	}

	/**  Move first period due to ascending trend 
	 @param serie the time series 
	 @return  the time series with moved first period 
	*/
	public final ArrayList<Double> MoveIniByUpTrend(ArrayList<Double> serie)
	{
		if (serie.size() <= minLeadTimeObs)
		{
			return serie;
		}
		ArrayList<Double> drv2 = poly.DeriveSerie(serie, 2);
		ArrayList<Double> trunk = serie;
		for (int index : poly.ZerosPosNegIndexes(drv2))
		{
			if (serie.size() - index < minLeadTimeObs)
			{
				return trunk;
			}
			trunk = new ArrayList<Double>();
			for (int i = index;i < serie.size();i++)
			{
				trunk.add(serie.get(i));
			}
			if (GetFunctionTrend(trunk) == 0)
			{
				return trunk;
			}


		}
		//if there are no more candidates
		int minIndex = 0;
		int maxIndex = trunk.size() - minLeadTimeObs;
		ArrayList<Double> tr = Trunk(trunk, maxIndex);

		//if it can be stationarized
		minIndex = GetMinIndex(trunk, 0, trunk.size() - 1, true, 0);
		trunk = Trunk(trunk, minIndex);
		return trunk;
	}

	/** Get index of minimmal data 
	 @param serie the time series 
	 @param minIndex start of interval 
	 @param maxIndex end of interval 
	 @param ascendant if the time series is ascendant 
	 @param it number of iterations 
	 @return  the index of minimmal data 
	*/
	public final int GetMinIndex(ArrayList<Double> serie, int minIndex, int maxIndex, boolean ascendant, int it)
	{
		if (maxIndex - minIndex <= 1 || it > 30)
		{
			return minIndex;
		}
		int res = ascendant? 1 : -1;
		int pivot = (minIndex + maxIndex) / 2;
		ArrayList<Double> trunk = Trunk(serie, pivot);
		if (GetFunctionTrend(trunk) == res)
		{
			return GetMinIndex(serie, pivot, maxIndex, ascendant, ++it);
		}
		else
		{
			return GetMinIndex(serie, minIndex, pivot, ascendant, ++it);
		}
	}

	//endregion

	//region Variance Trend

	/**  Calculate trend of the function 
	 @param serie the time series 
	 @return  value of trend 
	*/
	public final int GetFunctionTrend(ArrayList<Double> serie)
	{
		Result res = poly.WLSRegression(serie);
		double wlsTrend = res.c1; //Console.WriteLine("Trend: " + wlsTrend);
		if (wlsTrend < -trendDownThreshold)
		{
			return -1;
		}
		else if (wlsTrend > trendUpThreshold)
		{
			return 1;
		}
		else
		{
			return 0;
		}
	}

	/**  Calculate if the time series is stationary 
	 @param serie the time series 
	 @return  if it is stationary 
	*/
	public final boolean IsStationary(ArrayList<Double> serie)
	{
		Result res = poly.WLSRegression(serie);
		double wlsTrend = res.c1;
		return (wlsTrend >= -trendDownThreshold || wlsTrend <= trendUpThreshold);
	}

	/**  Calculate if the time series is homoscedastic 
	 @param serie the time series 
	 @return  if it is homoscedastic 
	*/
	public final boolean IsHomoscedastic(ArrayList<Double> serie)
	{
		resids = GetResiduals(serie);
		Result resResid = poly.WLSRegression(resids);
		double wlsResidTrend = resResid.c1;
		return (wlsResidTrend >= -trendDownThreshold || wlsResidTrend <= trendUpThreshold);
	}

	/**  Move the first period due to variance 
	 @param serie the time series 
	 @return  the time series with the first period moved 
	*/
	public final ArrayList<Double> MoveIniByVariance(ArrayList<Double> serie)
	{
		if (resids == null)
		{
			resids = GetResiduals(serie);
		}
		ArrayList<Double> drv2 = poly.DeriveSerie(resids, 2);
		ArrayList<Double> trunk = serie;
		for (int index : poly.ZerosNegPosIndexes(drv2))
		{
			if (serie.size() - index < minLeadTimeObs)
			{
				return trunk;
			}
			trunk = new ArrayList<Double>();
			for (int i = index;i < serie.size();i++)
			{
				trunk.add(serie.get(i));
			}
			if (GetFunctionTrend(trunk) == 0)
			{
				return trunk;
			}
		}
		return trunk;
	}

	private ArrayList<Double> GetResiduals(ArrayList<Double> serie)
	{
		Result res = poly.WLSRegression(serie);
		double wlsTrend = res.c1;
		ArrayList<Double> coeffs = new ArrayList<Double>();
		coeffs.add(res.c0);
		coeffs.add(res.c1);
		ArrayList<Double> x = new ArrayList<Double>();
		for (int i = 0;i < serie.size();i++)
		{
			x.add((double)i);
		}
		ArrayList<Double> wlsReg = poly.GetRegValues(x, coeffs);
		ArrayList<Double> resids = func.Substract(serie, wlsReg, true, false);
		return resids;
	}

	//endregion

	//region Filter Outliers

	/**  Filter outliers from the time series  
	 @param serie the time series 
	 @param outlierThreshold the outlier threshold 
	 @param minValueFiltered min value filtered, as a result 
	 @return  the filtered time series 
	*/
	public final ArrayList<Double> FilterOutliers(ArrayList<Double> serie, double outlierThreshold, DotNetHelpers.RefObject<Double> minValueFiltered)
	{
		Result res = poly.WLSRegression(serie);
		ArrayList<Double> coeffs = new ArrayList<Double>();
		coeffs.add(res.c0);
		coeffs.add(res.c1);
		ArrayList<Double> X = new ArrayList<Double>();
		for (int i = 0;i < serie.size();i++)
		{
			X.add((double)i);
		}
		ArrayList<Double> wlsReg = poly.GetRegValues(X, coeffs);

		ArrayList<Double> filtered = func.FilterOutliers(serie, wlsReg, outlierThreshold, minValueFiltered);
		return filtered;
	}

	//endregion

	//endregion

	//region Global Methods

	/**  Calculate first index to forecast 
	 @param serie the time series 
	 @return  the calculated first index 
	*/
	public final int GetFirstIndexToFcst(ArrayList<Double> serie)
	{
		boolean stationary = IsStationary(serie);
		boolean homoscedastic = IsHomoscedastic(serie);
		if (stationary && homoscedastic)
		{
			return 0;
		}
		else if (stationary && !homoscedastic)
		{
			return 0;
		}
		else
		{
			return 0;
		}
	}

	//endregion

	//region Private Methods

	private void Bisection(ArrayList<Double> serie, ArrayList<Double> initial, ArrayList<Double> final_Renamed)
	{
		int pivot = serie.size() / 2;
		for (int i = 0;i < pivot;i++)
		{
			initial.add(serie.get(i));
		}
		for (int i = pivot;i < serie.size();i++)
		{
			final_Renamed.add(serie.get(i));
		}
	}

	private TrendVar GetTrendVar(ArrayList<Double> serie)
	{
		Result res = poly.WLSRegression(serie);
		double wlsTrend = res.c1;
		double varCoeff = stat.VarCoeff(serie);
		if (wlsTrend < -trendDownThreshold)
		{
			if (varCoeff < varCoeffThreshold)
			{
				return TrendVar.TrendDesc_VarLow;
			}
			else
			{
				return TrendVar.TrendDesc_VarHigh;
			}
		}
		else if (wlsTrend > -trendDownThreshold)
		{
			if (varCoeff < varCoeffThreshold)
			{
				return TrendVar.TrendAsc_VarLow;
			}
			else
			{
				return TrendVar.TrendAsc_VarHigh;
			}
		}
		else
		{
			if (varCoeff < varCoeffThreshold)
			{
				return TrendVar.TrendEst_VarLow;
			}
			else
			{
				return TrendVar.TrendEst_VarHigh;
			}
		}
	}

	//endregion

	//region Enums

	/**  Result of trend-variance evaluation 
	*/
	public enum TrendVar
	{
		/**  trend ascendant, variance low 
		*/
		TrendAsc_VarLow,
		/**  trend stationary, variance low 
		*/
		TrendEst_VarLow,
		/**  trend descendant, variance low 
		*/
		TrendDesc_VarLow,
		/**  trend ascendant, variance high 
		*/
		TrendAsc_VarHigh,
		/**  trend stationary, variance high 
		*/
		TrendEst_VarHigh,
		/**  trend ascendant, variance high 
		*/
		TrendDesc_VarHigh;

		public int getValue()
		{
			return this.ordinal();
		}

		public static TrendVar forValue(int value)
		{
			return values()[value];
		}
	}

	//endregion

}