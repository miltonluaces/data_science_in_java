package Clustering;

import NumCalc.*;
import java.util.*;

/** 
 Descripci�n breve de Density.
*/
public class Density
{
	//region Fields 

	private ArrayList datos;
	private Functions func;

	//endregion

	//region Constructor

	/**  Constructor 
	 @param vectorDatos data matrix 
	*/
	public Density(double[][] vectorDatos)
	{
		func = new Functions();
		datos = new ArrayList();
		for (int i = 0;i < vectorDatos.length;i++)
		{
			datos.add(new Dato(vectorDatos[i]));
		}
		for (Dato dato : datos)
		{
			for (Dato otroDato : datos)
			{
				if (otroDato != dato)
				{
					dato.setDensidad(dato.getDensidad() + 1 / func.DistEuclid(dato.getValores(), otroDato.getValores()));
				}
			}
			dato.setDensidad(dato.getDensidad() / vectorDatos.length);
		}
	}

	//endregion

	//region Log Generation 

	/**  Generate log 
	*/
	public final void ClusteringLog()
	{
		try
		{
//C# TO JAVA CONVERTER WARNING: The java.io.OutputStreamWriter constructor does not accept all the arguments passed to the System.IO.StreamWriter constructor:
//ORIGINAL LINE: StreamWriter so = new StreamWriter("..\\Density.txt",false);
			java.io.OutputStreamWriter so = new java.io.OutputStreamWriter("..\\Density.txt");
			so.write("DENSITY: " + datos.size() + " Datos con " + ((Dato)datos.get(0)).getDim() + " elementos" + System.lineSeparator());
			so.write("--------------------------------------------------------------------------------" + System.lineSeparator());
			for (Dato dato : datos)
			{
				so.write(String.valueOf(dato) + System.lineSeparator());
			}
			so.close();
		}
		catch (java.lang.Exception e)
		{
			System.Diagnostics.Debug.WriteLine("Error");
		}
	}

	//endregion
}