package NumCalc;

import java.util.*;

/**  Class for assessing any forecast algorithms 
*/
public class FcstAssess
{

	//region Fields

	private String id;
	private FcstType type = FcstType.values()[0];
	private ArrayList<Double> fcst;
	private double dmd;

	private double stockLoss;
	private double stockLossRel;
	private double usDmd;
	private double usDmdRel;
	private double stockExcess;
	private double stockExcessRel;

	//endregion

	//region Constructor

	/**  Constructor 
	 @param id id of the assessment 
	 @param fcstType type of the forecast 
	*/
	public FcstAssess(String id, FcstType fcstType)
	{
		this.id = id;
		this.type = fcstType;
	}

	/**  Constructor 
	 @param id id of the assessment 
	 @param fcstType type of the forecast 
	 @param serie a list that contains the time series 
	 @param fcst a list htat contains the forecast 
	*/
	public FcstAssess(String id, FcstType fcstType, ArrayList<Double> serie, ArrayList<Double> fcst)
	{
		this.id = id;
		this.type = fcstType;
		this.fcst = fcst;
		Calculate(serie, fcst);
	}

	//endregion

	//region Properties

	/**  Id of the assess 
	*/
	public final String getId()
	{
		return id;
	}

	/**   Forecast type 
	*/
	public final FcstType getType()
	{
		return type;
	}

	/**   Stock loss metric 
	*/
	public final double getStockLoss()
	{
		return stockLoss;
	}
	public final void setStockLoss(double value)
	{
		stockLoss = value;
	}

	/**   Relative stock loss metric 
	*/
	public final double getStockLossRel()
	{
		return stockLossRel;
	}
	public final void setStockLossRel(double value)
	{
		stockLossRel = value;
	}

	/**   Unserved demand metric 
	*/
	public final double getUsDmd()
	{
		return usDmd;
	}
	public final void setUsDmd(double value)
	{
		usDmd = value;
	}

	/**  Relative unserved demand metric  
	*/
	public final double getUsDmdRel()
	{
		return usDmdRel;
	}
	public final void setUsDmdRel(double value)
	{
		usDmdRel = value;
	}

	/**  Excess in stock (for obsolescence) 
	*/
	public final double getStockExcess()
	{
		return stockExcess;
	}
	public final void setStockExcess(double value)
	{
		stockExcess = value;
	}

	/**  Relative stock excess  metric (for obsolescence) 
	*/
	public final double getStockExcessRel()
	{
		return stockExcessRel;
	}
	public final void setStockExcessRel(double value)
	{
		stockExcessRel = value;
	}

	/**  Demand 
	*/
	public final double getDmd()
	{
		return dmd;
	}
	public final void setDmd(double value)
	{
		dmd = value;
	}

	/**  Forecast 
	*/
	public final ArrayList<Double> getFcst()
	{
		return fcst;
	}

	//endregion

	//region Private Methods

	private void Calculate(ArrayList<Double> serie, ArrayList<Double> fcst)
	{
		dmd = CalcDmd(serie);
		stockLoss = CalcStockLoss(serie, fcst);
		stockLossRel = getStockLoss() / (double)serie.size();
		usDmd = CalcUsDmd(serie, fcst);
		if (dmd == 0.0)
		{
			usDmdRel = 0.0;
		}
		else
		{
			usDmdRel = usDmd / dmd;
		}
		stockExcess = CalcStockExcess(serie, fcst);
		if (dmd == 0.0)
		{
			stockExcessRel = 0.0;
		}
		else
		{
			stockExcessRel = getStockExcess() / dmd;
		}
	}

	private double CalcStockLoss(ArrayList<Double> serie, ArrayList<Double> fcst)
	{
		double sl = 0.0;
		for (int i = 0;i < serie.size();i++)
		{
			if (fcst.get(i).compareTo(serie.get(i)) < 0)
			{
				sl += 1.0;
			}
		}
		return sl;
	}

	private double CalcUsDmd(ArrayList<Double> serie, ArrayList<Double> fcst)
	{
		double usd = 0.0;
		for (int i = 0;i < serie.size();i++)
		{
			if (fcst.get(i).compareTo(serie.get(i)) < 0)
			{
				usd += serie.get(i) - fcst.get(i);
			}
		}
		return usd;
	}

	private double CalcStockExcess(ArrayList<Double> serie, ArrayList<Double> fcst)
	{
		double exc = 0.0;
		for (int i = 0;i < serie.size();i++)
		{
			if (fcst.get(i).compareTo(serie.get(i)) > 0)
			{
				exc += fcst.get(i) - serie.get(i);
			}
		}
		return exc;
	}

	private double CalcDmd(ArrayList<Double> fcst)
	{
		double dmd = 0.0;
		for (int i = 0;i < fcst.size();i++)
		{
			dmd += fcst.get(i);
		}
		return dmd;
	}

	//endregion

	//region enums

	/**  Types of forecast 
	*/
	public enum FcstType
	{
		/**  Forecast naïve (last value) 
		*/
		Naive,
		/**  Bayesian forecast 
		*/
		Bayes,
		/**  Leadtime forecast (by resampling, density composition, clustering; etc) 
		*/
		LtFcst,
		/**  Wavelet Neural Networks forecast (by wavelet frequency descomposition and additive neural network forecast 
		*/
		WnnFcst;

		public int getValue()
		{
			return this.ordinal();
		}

		public static FcstType forValue(int value)
		{
			return values()[value];
		}
	}

	//endregion

}