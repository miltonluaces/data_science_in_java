package NeuralNetworks;

import java.util.*;
import java.io.*;

/**  Class for reading txt files 
*/
public class ArchivoLectura
{

	//region Fields 

	private String nombreArchivo;
	private int cantFilas;
	private int maxColumnas;
	private char[] separador;

	//endregion

	//region Constructor

	/**  Constructor 
	 @param unN name 
	 @param unC max columns 
	 @param unS separation char 
	*/
	public ArchivoLectura(String unN, int unC, char unS)
	{
		nombreArchivo = unN;
		maxColumnas = unC;
		separador = new char[1];
		separador[0] = unS;
	}

	//endregion

	//region Properties

	/**  Number of rows  
	*/
	public final int getCantFilas()
	{
		return cantFilas;
	}

	//endregion

	//region Public Methods

	/**  Generate list of data from file 
	 @param cols number of columns 
	 @return  the list of data 
	*/
	public final ArrayList generarVector(int cols)
	{
		cantFilas = 0;
		ArrayList filas = new ArrayList();
		double[] columnas = new double[cols];
		String linea = "";
		String[] cifras;
		InputStreamReader arch = new java.io.InputStreamReader(nombreArchivo);
		linea = arch.ReadLine();
		double cifra;
		while (linea != null)
		{
			Debug.WriteLine("Linea : " + linea);
			columnas = new double[cols];
			cifras = linea.split(java.util.regex.Pattern.quote(separador.toString()), -1);
			cifra = Double.parseDouble(cifras[0]);
			columnas[0] = normalizar(cifra);
			for (int col = 1;col < cols;col++)
			{
				cifra = Double.parseDouble(cifras[col]);
				columnas[col] = (cifra == 0)? 0 : 1;
			}
			filas.add(columnas);
			cantFilas++;
			linea = arch.ReadLine();
		}
		return filas;
	}

	//endregion 

	//region Private Methods

	private double normalizar(double num)
	{
		return num / 6950.00; //revisar
	}

	//endregion

}