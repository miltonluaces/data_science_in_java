package Clustering;

import NumCalc.*;
import ClassicalStat.*;
import java.util.*;

//region Imports


//endregion

/**  Container class for SOM Clustering 
*/
public class Neurona
{

	//region Fields

	private int x;
	private int y;
	private int dim;
	private double[] w;
	private ArrayList datos;
	private double error;
	private Functions func;

	//endregion

	//region Constructor

	/**  Constructor 
	 @param dim number of dimensions 
	 @param x projection x axis width 
	 @param y projection y axis width 
	*/
	public Neurona(int dim, int x, int y)
	{
		this.func = new Functions();
		this.dim = dim;
		this.x = x;
		this.y = y;
		w = new double[dim];
		RandomGen rand = new RandomGen();
		rand.Reset();
		for (int i = 0;i < dim;i++)
		{
			w[i] = rand.NextDouble();
		}
	}

	//endregion

	//region Properties

	/**  projection x axis width 
	*/
	public final int getX()
	{
		return x;
	}

	/**  projection y axis width 
	*/
	public final int getY()
	{
		return y;
	}

	/**  number of dimensions 
	*/
	public final int getDim()
	{
		return dim;
	}
	public final void setDim(int value)
	{
		dim = value;
	}

	/**   arraylist of data 
	*/
	public final ArrayList getDatos()
	{
		return datos;
	}

	/**  global error 
	*/
	public final double getError()
	{
		return error;
	}
	public final void setError(double value)
	{
		error = value;
	}

	//endregion

	//region Setters and Getters

	/**  Increment weight 
	 @param i index 
	 @param peso weight 
	*/
	public final void IncW(int i, double peso)
	{
		w[i] += peso;
	}

	/**  Get a weight from an index 
	 @param i the index 
	 @return  the weight 
	*/
	public final double GetW(int i)
	{
		return w[i];
	}

	//endregion

	//region Public Methods

	/**  Euclidian distance 
	 @param valores array of values 
	 @return  the distance value 
	*/
	public final double Distancia(double[] valores)
	{
		return func.DistEuclid(w, valores);
	}

	/**  Add new datum 
	 @param dato the datum to add 
	*/
	public final void AddDato(Dato dato)
	{
		datos.add(dato);
	}

	/**  Clear the neuron from data 
	*/
	public final void Reset()
	{
		datos = new ArrayList();
	}

	//endregion

	//region ToString Override

	/**  ToString override  
	 @return  the formatted string 
	*/
	@Override
	public String toString()
	{
		String ret = "[" + x + "," + y + "] : \t";
		for (Dato dato : datos)
		{
			ret += dato + " ";
		}
		return ret;
	}

	//endregion
}