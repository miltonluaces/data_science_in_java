package SignalProcess;

import java.util.*;

/**  Abstract class for any signal processing method 
*/
public abstract class SignalProc
{

	//region Fields

	/**   Original time series 
	*/
	protected ArrayList<Double> origData;

	/**   Current time series 
	*/
	protected ArrayList<Double> data;

	/**   Trend time series 
	*/
	protected ArrayList<Double> trend;

	/**  proportion of non-zeros 
	*/
	protected double loadFactor;


	//endregion

	//region Properties

	/**  Original time series 
	*/
	public final ArrayList<Double> getOrigData()
	{
		return origData;
	}

	/**  Current time series 
	*/
	public final ArrayList<Double> getData()
	{
		return data;
	}

	/**  Time series trend at each point 
	*/
	public final ArrayList<Double> getTrend()
	{
		return trend;
	}

	/**  load factor of grouped non-zero values on total grouped values 
	*/
	public final double getLoadFactor()
	{
		return loadFactor;
	}

	//endregion

	//region Public Methods

	/**  Group time series data in groups of n values 
	 @param origData original time series 
	 @param firstIndex index of first period for calculation 
	 @param group grouping index 
	 @param continuous if applies continuous grouping or not 
	*/
	public final void Group(ArrayList<Double> origData, int firstIndex, int group, boolean continuous)
	{
		if (continuous)
		{
			ContinuousGroup(origData, firstIndex, group);
		}
		else
		{
			SimpleGroup(origData, firstIndex, group);
		}
	}

	public final void SimpleGroup(ArrayList<Double> origData, int firstIndex, int group)
	{
		if (group == 1)
		{
			this.trend = new ArrayList<Double>(origData);
		this.data = origData;
			return;
		}

		this.origData = origData;
		this.data = new ArrayList<Double>();
		int day = 0;
		double totPeriodo = 0.0;
		double val;
		int nonZero = 0;
		for (int i = origData.size() - 1;i >= firstIndex;i--)
		{
			val = origData.get(i);
			totPeriodo += val;
			day++;
			if (day == group)
			{
				data.add(0,totPeriodo);
				day = 0;
				if (totPeriodo > 0)
				{
					nonZero++;
				}
				totPeriodo = 0.0;
			}
		}
		if (day > 0)
		{
			data.add(0,totPeriodo);
			if (totPeriodo > 0)
			{
				nonZero++;
			}
		}
		loadFactor = (double)nonZero / (double)data.size();
		this.trend = new ArrayList<Double>(data);
	}

	public final void ContinuousGroup(ArrayList<Double> origData, int firstIndex, int group)
	{
		this.data = new ArrayList<Double>();
		if (origData == null || origData.isEmpty() || group <= 0)
		{
			return;
		}
		int nonZero = 0;
		for (int i = firstIndex; i < origData.size() - group; i++)
		{
			double total = 0;

			int l = 0;
			int j = i;
			while (l < group)
			{
				if (j >= origData.size())
				{
					return;
				}
				if (origData.get(j).compareTo(0) >= 0)
				{
					total += origData.get(j);
				}
				l++;
				j++;
			}
			data.add(total);
			if (total > 0)
			{
				nonZero++;
			}
		}
		loadFactor = (double)nonZero / (double)data.size();
		this.trend = new ArrayList<Double>(data);
	}


	//endregion

}