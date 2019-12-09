package NumCalc;

import java.util.*;

/**  Obsolete class that represents a distribution (Histogram is recommended instead of this class) 
*/
public class Distribution
{

	//region Fields

	private TreeMap<Integer, Column> columns;
	private int count;
	private int exp;
	private HashMap<Integer, Double> percentiles;

	//endregion

	//region Constructor

	/**  Constructor 
	*/
	public Distribution()
	{
		columns = new TreeMap<Integer, Column>();
		percentiles = new HashMap<Integer, Double>();
	}

	//endregion

	//region Properties

	/**  Dictionary that contains the column objects 
	*/
	public final TreeMap<Integer, Column> getColumns()
	{
		return columns;
	}

	/**  number of data 
	*/
	public final int getCount()
	{
		return count;
	}

	/**  Exponent (Power of 10) Positive or negative (number of zeros) 
	*/
	public final int getExp()
	{
		return exp;
	}

	//endregion

	//region Public Methods

	/**  Load data from a list of values 
	 @param values the data values 
	 @param exp exponent 
	*/
	public final void LoadValues(ArrayList<Double> values, int exp)
	{
		this.count = values.size();
		this.exp = exp;
		this.getColumns().clear();
		double discVal = 0.0;

		for (double val : values)
		{
			if (exp < 0)
			{
				discVal = Math.round(val * Math.pow(10, -exp)) / Math.pow(10, -exp);
			}
			else
			{
				discVal = Math.round((val / Math.pow(10,exp)) * Math.pow(10, 0)) / Math.pow(10, 0);
			}
			if (!Contains(discVal))
			{
				Add(discVal);
			}
			IncFrec(discVal);
		}
		int acum = 0;
		for (Column col : columns.values())
		{
			acum = acum + col.freq;
			col.acum = acum;
		}
		double totAcum = (double)acum;
		for (Column col : columns.values())
		{
			col.prob = (double)col.acum / totAcum;
		}
	}

	/**  Load data from a list that represents a distribution (indexes = values, values = frequencies) 
	 @param meanFrec the distribution list 
	 @param exp exponent 
	*/
	public final void LoadDist(ArrayList<Double> meanFrec, int exp)
	{
		this.count = meanFrec.size();
		this.exp = exp;
		this.getColumns().clear();
		for (int i = 0;i < meanFrec.size();i++)
		{
			Add((double)i);
			SetFrec((double)i, meanFrec.get(i));
		}
		int acum = 0;
		for (Column col : columns.values())
		{
			acum = acum + col.freq;
			col.acum = acum;
		}
		double totAcum = (double)acum;
		for (Column col : columns.values())
		{
			col.prob = (double)col.acum / totAcum;
		}
	}

	/**  Calculate percentile for any probability 
	 @param prob the probability 
	 @return  the value of the percentile 
	*/
	public final double CalculatePercentile(double prob)
	{
		double val = 0;
		Column col;
		for (int key : columns.keySet())
		{
			col = columns.get(key);
			if (col.prob > prob)
			{
				return val;
			}
			val = col.val;
		}
		return val;
	}

	/**  Get a calculated percentile of any probability 
	 @param p the probability 
	 @return  the calculated value for the percentile 
	*/
	public final double GetPercentile(double p)
	{
		return percentiles.get(GetKey(p));
	}

	/**  Get the probability for any value  
	 @param x the value 
	 @return  the calculated probability 
	*/
	public final double GetProbability(double x)
	{
		Column prevCol = null;
		for (int key : columns.keySet())
		{
			if (columns.get(key).val > x)
			{
				if (prevCol == null || columns.get(key).val - x < columns.get(key).val - prevCol.val)
				{
					return columns.get(key).prob;
				}
				else
				{
					return prevCol.prob;
				}
			}
			prevCol = columns.get(key);
		}
		return 1.00;
	}

	/**  Get a probability of an interval 
	 @param lower lower limit of the interval 
	 @param upper upper limit of the interval 
	 @return  the calculated probability 
	*/
	public final double GetProbability(double lower, double upper)
	{
		double lowProb = GetProbability(lower);
		double uppProb = GetProbability(upper);
		return upper - lower;
	}

	/**  Set a percentile 
	 @param p the percentage 
	 @param val the value to set 
	*/
	public final void SetPercentil(double p, double val)
	{
		percentiles.put(GetKey(p), val);
	}

	/**  Set all percentiles 
	*/
	public final void SetPercentiles()
	{
		percentiles.clear();
		double p, x;
		double xAnt = -Double.MAX_VALUE;
		for (int i = 800;i < 1000;i++)
		{
			p = (double)i / 10.0;
			x = CalculatePercentile(p * 0.01);
			if (x < xAnt)
			{
				x = xAnt;
			}
			percentiles.put(GetKey(p), x);
			xAnt = x;
		}
	}

	//endregion

	//region Private Methods

	private Column GetColumn(double val)
	{
		return columns.get(GetKey(val));
	}

	private boolean Contains(double val)
	{
		return columns.containsKey(GetKey(val));
	}

	private void Add(double val)
	{
		columns.put(GetKey(val), new Column(val));
	}

	private void IncFrec(double val)
	{
		columns.get(GetKey(val)).IncFrec();
	}

	private void SetFrec(double val, double freq)
	{
		columns.get(GetKey(val)).freq = (int)java.lang.Math.round(Math.round(freq));
	}

	private int GetKey(double val)
	{
		return (int)java.lang.Math.round(Math.round(val * Math.pow(10, 1)) / Math.pow(10, 1) * 10);
	}

	/**  Get correspondent value of a particular key  
	 @param key the int key 
	 @return  the corresponding value 
	*/
	public final double GetValue(int key)
	{
		return (double)key / 10.0;
	}

	private void AddPercentile(double perc, double val)
	{
		percentiles.put(GetKey(perc), val);
	}


	//endregion

	//region Inner class Column

	/**  Column of the histogram  
	*/
	public static class Column
	{

		/**  value 
		*/
		public double val;
		/**  frequency 
		*/
		public int freq;
		/**  acumulated frequency 
		*/
		public int acum;
		/**  probability 
		*/
		public double prob;

		/**  Constructor 
		 @param val value 
		*/
		public Column(double val)
		{
			this.val = val;
			this.freq = 0;
		}

		/**  Increment frequency 
		*/
		public final void IncFrec()
		{
			this.freq = this.freq + 1;
		}
	}

	//endregion

}
