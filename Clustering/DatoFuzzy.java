package Clustering;

import java.util.*;

/**  A wrapper class for a fuzzy datum 
*/
public class DatoFuzzy extends Dato
{

	//region Fields

	private int nClusters;
	private double[] M;

	//endregion

	//region Constructors

	/**  Constructor 
	 @param dim number of dimensions 
	*/
	public DatoFuzzy(int dim)
	{
		super(dim);
	}

	/**  Constructor 
	 @param valores array of values to load 
	*/
	public DatoFuzzy(double[] valores)
	{
		super(valores);
	}

	/**  Constructor 
	 @param dim number of dimensions 
	 @param nClusters number of clusters 
	*/
	public DatoFuzzy(int dim, int nClusters)
	{
		super(dim);
		this.nClusters = nClusters;
		M = new double[nClusters];
		double totM = 0f;
		for (int i = 0;i < nClusters;i++)
		{
			M[i] = (new Random()).nextInt();
			totM += M[i];
		}
		for (int i = 0;i < nClusters;i++)
		{
			M[i] = M[i] / totM; // Condici�n: � M[i] = 1
		}
	}

	/**  Constructor 
	 @param valores array of values to load 
	 @param nClusters number of clusters 
	*/
	public DatoFuzzy(double[] valores, int nClusters)
	{
		super(valores);
		this.nClusters = nClusters;
		M = new double[nClusters];
		double totM = 0f;
		for (int i = 0;i < nClusters;i++)
		{
			M[i] = (new Random()).nextInt();
			totM += M[i];
		}
		for (int i = 0;i < nClusters;i++)
		{
			M[i] = M[i] / totM; // Condici�n: � M[i] = 1
		}
	}

	//endregion

	//region Setters and Getters

	/**  Set a fuzzy value in membership matrix 
	 @param valor the fuzzy value 
	 @param pos position 
	*/
	public final void SetM(double valor, int pos)
	{
		M[pos] = valor;
	}

	/**  Get a fuzzy value from membership matrix 
	 @param pos position 
	 @return  the fuzzy value 
	*/
	public final double GetM(int pos)
	{
		return M[pos];
	}

	//endregion

	//region To String Override

	/**  ToString override  
	 @return  the formatted string 
	*/
	@Override
	public String toString()
	{
		String str = super.toString();
		if (M != null)
		{
			str += " - M[";
			for (int i = 0;i < nClusters;i++)
			{
				str += M[i] + " ";
			}
			str += "]";
		}
		return str;
	}

	//endregion

}