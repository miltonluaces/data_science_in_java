package NumCalc;


import java.util.*;

/**  Frecuencies of a set of values 
*/
public class Frequencies
{

	//region Fields

	private ArrayList<Double> frec;
	private int nValues;
	private double min;
	private double max;
	private double sum;
	private double sum2;

	//endregion

	//region Constructor

	/**  Constructor 
	*/
	public Frequencies()
	{
		frec = new ArrayList<Double>();
		Clear();
	}

	//endregion

	//region Load Data

	/**  Clears all frequencies 
	*/
	public final void Clear()
	{
		frec.clear();
		nValues = 0;
		min = Double.MAX_VALUE;
		max = -1;
		sum = 0.0;
		sum2 = 0.0;
	}

	/**  Adds a value 
	 @param valor the value 
	*/
	public final void Add(double valor)
	{
		int val = (int)valor;
		if (val < min)
		{
			min = val;
		}
		if (val >= frec.size() - 1)
		{
			for (int i = frec.size();i <= val + 1;i++)
			{
				frec.add(0.0);
			}
			max = val;
		}
		frec.set(val, frec.get(val) + 1);
		nValues++;
		sum += val;
		sum2 += Math.pow(val, 2);
	}

	/**  Adds an absolute value  
	 @param valor the value 
	*/
	public final void AddAbs(double valor)
	{
		Add(Math.abs(valor));
	}

	/**  Removes a value 
	 @param valor the value 
	*/
	public final void Delete(double valor)
	{
		int val = (int)valor;
		if (frec.get(val).equals(0))
		{
			throw new RuntimeException("Cannot_delete_values_with_zero_frequency");
		}

		frec.set(val, frec.get(val) - 1);
		nValues--;
		sum -= val;
		sum2 -= Math.pow(val, 2);

		if (val == min && frec.get(val).equals(0))
		{
			int newMin = val + 1;
			while (frec.get(newMin).equals(0))
			{
				newMin++;
			}
			min = (double)newMin;
		}
		if (val == max && frec.get(val).equals(0))
		{
			int newMax = val - 1;
			while (newMax == 0)
			{
				frec.remove(newMax);
				newMax--;
			}
			max = (double)newMax;
		}
	}

	/**  Adds n values 
	 @param valores an array of value 
	*/
	public final void AddRange(double[] valores)
	{
		for (int val : valores)
		{
			Add(val);
		}
	}

	/**  Adds a value as double 
	 @param val the int value 
	*/
	public final void Add(int val)
	{
		Add((double)val);
	}

	/**  Adds n values in a IList of ints 
	 @param valores IList of ints 
	*/
	public final void AddRange(List<Integer> valores)
	{
		for (int val : valores)
		{
			Add(val);
		}
	}

	/**  Adds n values ina a IList of doubles 
	 @param valores IList of doubles 
	*/
	public final void AddRange(List<Double> valores)
	{
		for (int val : valores)
		{
			Add(val);
		}
	}

	/**  Add n values in a IList of doubles in absolute value 
	 @param valores IList of doubles 
	*/
	public final void AddRangeAbs(List<Double> valores)
	{
		for (int val : valores)
		{
			AddAbs(val);
		}
	}

	/**  Removes n values in a IList of doubles 
	 @param valores IList of doubles 
	*/
	public final void DeleteRange(List<Double> valores)
	{
		for (int val : valores)
		{
			Delete(val);
		}
	}


	/**  Normalized data load 
	 @param valores List of values 
	 @param cotaSup superior cota 
	*/
	public final void LoadNormalizedData(ArrayList<Double> valores, double cotaSup)
	{
		Clear();
		ArrayList<Double> valoresNorm = new ArrayList<Double>();
		Norma norma = new Norma(valores);
		norma.Normalize(Norma.NormType.minMax);
		if (norma.getMax() <= cotaSup)
		{
			valoresNorm = valores;
		}
		else
		{
			Hashtable clases = new Hashtable();
			Hashtable medias = new Hashtable();
			double key;
			for (int i = 0;i < valores.size();i++)
			{
				key = (int)(cotaSup * norma.GetValue(i, true));
				if (clases.get(key) == null)
				{
					clases.put(key, new ArrayList<Double>());
				}
				((ArrayList<Double>)clases.get(key)).add(valores.get(i));
			}
			for (double keY : clases.keySet())
			{
				double tot = 0.0;
				ArrayList<Double> values = (ArrayList<Double>)clases.get(keY);
				for (double val : values)
				{
					tot += val;
				}
				double prom = tot / values.size();
				medias.put(keY, cotaSup * ((prom - norma.getMin()) / norma.getRango()));
			}

			valoresNorm = new ArrayList<Double>();
			for (int i = 0;i < valores.size();i++)
			{
				key = (int)(cotaSup * norma.GetValue(i, true));
				valoresNorm.add((double)medias.get(key));
			}
		}
		AddRange(valoresNorm);
	}

	/**  Data load from a dictionary  
	 @param dict the dictionary 
	*/
	public final void LoadData(HashMap<Integer, Double> dict)
	{

		Clear();
		double val;
		for (int key : dict.keySet())
		{
			val = dict.get(key);
			if (key < min)
			{
				min = key;
			}
			nValues += (int)val;
			sum += key * val;
			sum2 += Math.pow(key, 2) * val;
			if (key >= frec.size() - 1)
			{
				for (int i = frec.size();i <= key + 1;i++)
				{
					frec.add(0.0);
				}
				max = key;
			}
			frec.set(key, val);
		}
	}

	/**  Returns a collection of values repeated any times acording to each frequency/ 
	 @return  the list of repeated values 
	*/
	public final ArrayList<Double> GetValuesXFrec()
	{
		ArrayList<Double> valXFrec = new ArrayList<Double>();
		for (int v = 0;v < this.frec.size();v++)
		{
			for (int f = 0;f < (int)frec.get(v);f++)
			{
				valXFrec.add(v);
			}
		}
		return valXFrec;
	}

	/**  Loads a dictionary of values and its not-null frequencies 
	 @param valXFrec values per frequency 
	*/
	public final void GetValuesXFrecDict(HashMap<Integer, Double> valXFrec)
	{
		for (int i = 0;i < this.frec.size();i++)
		{
			if (!frec.get(i).equals(0))
			{
				valXFrec.put(i, frec.get(i));
			}
		}
	}

	//endregion

	//region Properties

	/**  Minimum 
	*/
	public final double getMin()
	{
		return min;
	}

	/**  Maximum 
	*/
	public final double getMax()
	{
		return max;
	}

	/**  Range 
	*/
	public final double getRange()
	{
		return max - min;
	}

	/**  Number of values 
	*/
	public final int getNValues()
	{
		return nValues;
	}

	/**  Number of frequencies 
	*/
	public final int getNFrecs()
	{
		return frec.size();
	}

	/**  Frequencies (included zero in max + 1) 
	*/
	public final ArrayList<Double> getValFrecs()
	{
		return frec;
	}

	/**  Sum of values 
	*/
	public final double getSum()
	{
		return sum;
	}

	/**  Squared sum of values 
	*/
	public final double getSum2()
	{
		return sum2;
	}

	/**  Values mean 
	*/
	public final double getMean()
	{
		return getSum() / getNValues();
	}

	/**  Values variance 
	*/
	public final double getVar()
	{
		return (getNValues() * sum2 - Math.pow(sum, 2)) / (getNValues() * (getNValues() - 1.0));
	}

	/**  Standard deviation 
	*/
	public final double getStDev()
	{
		return Math.sqrt(getVar());
	}

	/**  Populates an arrayList of values 
	*/
	public final ArrayList getValues()
	{
		ArrayList values = new ArrayList();
		for (int i = 0;i < getNFrecs();i++)
		{
			values.add((double)i);
		}
		return values;
	}

	/**  Populates an arrayList of frequencies 
	*/
	public final ArrayList getFrecs()
	{
		ArrayList frecs = new ArrayList();
		for (int i = 0;i < getNFrecs();i++)
		{
			frecs.add((double)frec.get(i));
		}
		return frecs;
	}

	//endregion

	//region ToString Override

	/**  ToString override  
	 @return  the formatted string 
	*/
	@Override
	public String toString()
	{
		String str = Properties.Strings.Frecuencias;
		str += (Properties.Strings.Valor___frecuencia);
		for (int i = 0;i < getNFrecs();i++)
		{
			str += (i + "\t" + frec.get(i).toString("0.000")) + "\n";
		}
		return str;
	}

	/** 
	 To String for not null frequencies
	 
	 @return  formatted string 
	*/
	public final String ToStringNotNull()
	{
		String str = Properties.Strings.FrecuenciasNoNulas;
		str += (Properties.Strings.Valor___frecuencia);
		for (int i = 0;i < getNFrecs();i++)
		{
			if (!frec.get(i).equals(0))
			{
				str += (i + "\t" + frec.get(i).toString()) + "\n";
			}
		}
		return str;
	}

	//endregion

}