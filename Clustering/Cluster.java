package Clustering;

import NumCalc.*;
import java.util.*;

/**  Class for different kind of Clustering methods  
*/
public class Cluster implements Comparator, java.lang.Comparable
{

	//region Fields

	private Functions func;
	private double indice;
	private Dato centro;
	private ArrayList datos;
	private double dif;
	private double dist;

	//endregion

	//region Constructors

	/**  Constructor 
	*/
	public Cluster()
	{
		func = new Functions();
		centro = null;
		datos = new ArrayList();
		dif = 0.0;
	}

	/**  Constructor 
	 @param centro centre of the cluster 
	*/
	public Cluster(Dato centro)
	{
		func = new Functions();
		this.centro = centro;
		datos = new ArrayList();
		dif = 0.0;
	}

	//endregion

	//region Properties

	/**   index of the cluster 
	*/
	public final double getIndice()
	{
		return indice;
	}
	public final void setIndice(double value)
	{
		indice = value;
	}

	/**   centre of the cluster 
	*/
	public final Dato getCentro()
	{
		return centro;
	}
	public final void setCentro(Dato value)
	{
		centro = value;
	}

	/**   list of data in this cluster 
	*/
	public final ArrayList getDatos()
	{
		return datos;
	}

	/**   quantity of data 
	*/
	public final int getNDatos()
	{
		return datos.size();
	}

	/**   difference 
	*/
	public final double getDif()
	{
		return dif;
	}
	public final void setDif(double value)
	{
		dif = value;
	}

	/**   Discrimination 
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

	//region Setters and Getters

	/**  Get data by index 
	 @param i index 
	 @return  the requested datum 
	*/
	public final Dato GetDato(int i)
	{
		return (Dato)datos.get(i);
	}

	/**  Add a new datum 
	 @param unD the datum to add 
	*/
	public final void Add(Dato unD)
	{
		datos.add(unD);
	}

	/**   Add a range of data 
	 @param datos the data in an ArrayList 
	*/
	public final void Add(ArrayList datos)
	{
		for (Dato dato : datos)
		{
			Add(dato);
		}
	}

	//endregion

	//region Public Methods

	/**  Calculate cluster centre 
	 @return  the centre 
	*/
	public final Dato CalcularCentro()
	{
		int dim = centro.getDim();
		if (getNDatos() > 0)
		{ //          \� K            [df,K = cant datos]
														//          /_ df=1 U kd     [d,D = dims]
			double[] totDim = new double[dim]; //  Ce =  ���������������   [U kv = valor v del dato df]
			for (int d = 0; d < dim; d++)
			{
				for (int k = 0; k < getNDatos(); k++)
				{ //              K
					totDim[d] += ((Dato)datos.get(k))[d];
				}
			}

			for (int d = 0; d < dim; d++)
			{
				totDim[d] = totDim[d] / (double)getNDatos();
			}
			centro = new Dato(totDim);
			return centro;
		}
		else
		{
			return centro;
		}
	}

	/**  Mean distance of data to the cluster centre 
	 @return  the distance 
	*/
	public final double DistanciaMedia()
	{
		double totDist = 0f;
		for (int i = 0;i < getNDatos();i++)
		{
			totDist += func.DistEuclid(centro.getValores(), ((Dato)datos.get(i)).getValores());
		}
		if (getNDatos() == 0)
		{
			return 0.0;
		}
		return totDist / getNDatos();
	}

	/**  Clear cluster data values 
	*/
	public final void Reset()
	{
		datos = new ArrayList();
		dif = 0.0;
	}

	//endregion

	//region ToString Override

	/**  ToString override 
	 @return  the string to show 
	*/
	@Override
	public String toString()
	{
		String ret = "";
		ret += String.format(Properties.Strings.Centro_datos, centro, getNDatos());
		for (int i = 0;i < getNDatos();i++)
		{
			ret += datos.get(i) + " ";
		}
		//ret += "\nDistancia Media: "+ DistanciaMedia();
		return ret;
	}

	//endregion

	//region Implementacion de IComparable 

	public final int compareTo(Object obj)
	{
		return (int)(this.getIndice() * 1000000 - ((Cluster)obj).getIndice() * 1000000);
	}

	public final int compare(Object x, Object y)
	{
		if (((Cluster)x).getIndice() * 1000000 > ((Cluster)y).getIndice() * 1000000)
		{
			return -1;
		}
		else if (((Cluster)x).getIndice() * 1000000 < ((Cluster)y).getIndice() * 1000000)
		{
			return 1;
		}
		else
		{
			return 0;
		}
	}

	//endregion

	//region Sort by Property Inner Classes

	private static class sortDistAscHelper implements Comparator<Cluster>
	{
		public final int compare(Cluster c1, Cluster c2)
		{
			if (c1.getDist() < c2.getDist())
			{
				return 1;
			}
			if (c1.getDist() > c2.getDist())
			{
				return -1;
			}
			else
			{
				return 0;
			}
		}
	}


	private static class sortDistDescHelper implements Comparator<Cluster>
	{
		public final int compare(Cluster c1, Cluster c2)
		{
			if (c1.getDist() < c2.getDist())
			{
				return 1;
			}
			if (c1.getDist() > c2.getDist())
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

	public static Comparator<Cluster> sortDistAsc()
	{
		return (Comparator<Cluster>)new sortDistAscHelper();
	}
	public static Comparator<Cluster> sortDistDesc()
	{
		return (Comparator<Cluster>)new sortDistDescHelper();
	}

	//endregion

}