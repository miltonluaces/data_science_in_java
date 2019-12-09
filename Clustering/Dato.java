package Clustering;

import NumCalc.*;
import java.util.*;

/**  Container class representing a single datum (datum wrapper) 
*/
public class Dato implements java.lang.Comparable
{

	//region Fields

	private String codigo;
	private int dim;
	private int sortDim;
	private double[] valores;
	private double[] norm;
	private double densidad;
	private Cluster cluster;
	private double dist;

	//endregion

	//region Constructors

	/**  Constructor 
	 @param dim number of dimensions 
	*/
	public Dato(int dim)
	{
		this.dim = dim;
		this.valores = new double[dim];
		this.norm = new double[dim];
		for (int i = 0;i < dim;i++)
		{
			valores[i] = (new Random()).nextInt();
		}
	}

	/**  Constructor 
	 @param valores array of values to load 
	*/
	public Dato(double[] valores)
	{
		this.valores = valores;
		dim = valores.length;
		this.norm = new double[dim];
	}

	/**  Constructor 
	 @param valores array of values to load 
	 @param pesos array of weights to load 
	*/
	public Dato(double[] valores, double[] pesos)
	{
		this.valores = valores;
		dim = valores.length;
		this.norm = new double[dim];
	}

	/**  Constructor 
	 @param codigo datum code 
	 @param valores array of values to load 
	 @param pesos array of weights to load 
	*/
	public Dato(String codigo, double[] valores, double[] pesos)
	{
		this.codigo = codigo;
		this.valores = valores;
		dim = valores.length;
		this.norm = new double[dim];
	}

	/**  Constructor 
	 @param codigo datum code 
	 @param valores array of values to load 
	*/
	public Dato(String codigo, double[] valores)
	{
		this.codigo = codigo;
		this.valores = valores;
		dim = valores.length;
		this.norm = new double[dim];
	}

	//endregion

	//region Properties

	/**  Datum code 
	*/
	public final String getCodigo()
	{
		return codigo;
	}
	public final void setCodigo(String value)
	{
		codigo = value;
	}

	/**  Number of dimensions 
	*/
	public final int getDim()
	{
		return dim;
	}
	public final void setDim(int value)
	{
		dim = value;
	}

	/**  Number of dimension sorted 
	*/
	public final int getSortDim()
	{
		return sortDim;
	}
	public final void setSortDim(int value)
	{
		sortDim = value;
	}

	/**  Cluster from where the datum belongs 
	*/
	public final Cluster getCluster()
	{
		return cluster;
	}
	public final void setCluster(Cluster value)
	{
		cluster = value;
	}

	/**  Density of this datum 
	*/
	public final double getDensidad()
	{
		return densidad;
	}
	public final void setDensidad(double value)
	{
		densidad = value;
	}

	/**  Values of each dimension of the datum 
	*/
	public final double[] getValores()
	{
		return valores;
	}

	/**   Normalized values of each dimension of the datum 
	*/
	public final double[] getNorm()
	{
		return norm;
	}

	/**   Indexer of the datum 
	*/
	public final double getItem(int i)
	{
		return valores[i];
	}
	public final void setItem(int i, double value)
	{
		valores[i] = value;
	}

	/**  Dist for discrimination 
	*/
	public final double getDist()
	{
		return dist;
	}
	public final void setDist(double value)
	{
		dist = value;
	}

	//endregion

	//region To String Override

	/**  ToString override 
	 @return  the formatted string 
	*/
	@Override
	public String toString()
	{
		String str = ""; //"[";
		for (int i = 0; i < dim; i++)
		{
			str += (new Double(valores[i])).toString("////.00") + "\t";
		}
		//str += "] ";
		if (densidad != 0)
		{
			str += " (dens = " + (new Double(densidad)).toString("////.00") + ")";
		}
		return str;
	}

	//endregion

	//region Public Methods 

	/**  Clone a new datum from this (deep copy) 
	 @return  the cloned datum 
	*/
	public final Dato Clone()
	{
		Dato clon = new Dato(this.dim);
		clon.codigo = this.codigo;
		clon.sortDim = this.sortDim;
		clon.densidad = this.densidad;
		clon.valores = new double[dim];
		for (int i = 0;i < dim;i++)
		{
			clon.valores[i] = this.valores[i];
		}
		clon.setCluster(new Cluster());
		return clon;
	}

	/**  equal operator override 
	 @param dato1 first datum 
	 @param dato2 second datum 
	 @return  if they are equal 
	*/
	public static boolean OpEquality(Dato dato1, Dato dato2)
	{
		if (dato1.dim != dato2.dim)
		{
			return false;
		}
		for (int i = 0;i < dato1.dim;i++)
		{
			if (dato1.getValores()[i] != dato2.getValores()[i])
			{
				return false;
			}
		}
		return true;
	}

	/**  not equal operator override 
	 @param dato1 first datum 
	 @param dato2 second datum 
	 @return  if they are not equal 
	*/
	public static boolean OpInequality(Dato dato1, Dato dato2)
	{
		return !(dato1 == dato2);
	}

	/**  Get Datum code 
	 @return  the code 
	*/
	public final String GetCode()
	{
		String code = "";
		for (double val : getValores())
		{
			code += val;
		}
		return code;
	}

	//endregion

	//region IComparable Members

	public final int compareTo(Object obj)
	{
		return (int)(this.getValores()[sortDim] - ((Dato)obj).getValores()[sortDim]);
	}

	//endregion

	//region Sort by Property Inner Classes

	private static class sortDistAscHelper implements Comparator<Dato>
	{
		public final int compare(Dato d1, Dato d2)
		{
			if (d1.getDist() < d2.getDist())
			{
				return 1;
			}
			if (d1.getDist() > d2.getDist())
			{
				return -1;
			}
			else
			{
				return 0;
			}
		}
	}


	private static class sortDistDescHelper implements Comparator<Dato>
	{
		public final int compare(Dato d1, Dato d2)
		{
			if (d1.getDist() < d2.getDist())
			{
				return 1;
			}
			if (d1.getDist() > d2.getDist())
			{
				return -1;
			}
			else
			{
				return 0;
			}
		}
	}

	//endregion

	//region Sort by Property Static Methods

	public static Comparator<Dato> sortDistAsc()
	{
		return (Comparator<Dato>)new sortDistAscHelper();
	}
	public static Comparator<Dato> sortDistDesc()
	{
		return (Comparator<Dato>)new sortDistDescHelper();
	}

	//endregion

	//region Overrides

	/**  Equals override
	 @param obj allocation to compare 
	 @return  if objects are equal 
	*/
	@Override
	public boolean equals(Object obj)
	{
		return super.equals(obj);
	}

	/**  Override for hash function 
	 @return  hashcode 
	*/
	@Override
	public int hashCode()
	{
		return super.hashCode();
	}

	//endregion

}