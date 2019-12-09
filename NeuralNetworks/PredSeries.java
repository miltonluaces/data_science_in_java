package NeuralNetworks;

import java.util.*;

//region Imports


//endregion

/**  Class Predictability of time series (supervised learning) 
*/
public class PredSeries
{

	//region Fields

	private Red [] Red;
	private int [] Size;

	//endregion

	//region Singleton

	private static PredSeries instance;

	/**  The method that returns the one and only instance 
	 @param sizeDaily size of time series in case of dayly periodicity 
	 @param n trainning coefficent 
	 @param m moment factor 
	 @return  the one and only instance 
	*/
	public static PredSeries getInstance(int sizeDaily, double n, double m)
	{
		if (instance == null)
		{
			instance = new PredSeries(sizeDaily, n, m);
		}
		return instance;
	}

	//endregion

	//region Private Constructor

	private PredSeries(int sizeDaily, double n, double m)
	{
		Size = new int[4];
		Size[1] = 120; //Monthly
		Size[2] = 265; //Weekly
		Size[3] = sizeDaily;

		Red = new Red[4];
		for (int i = 1;i <= 3;i++)
		{
			Red[i] = new Red(i,Size[i] + 1,(int)(Size[i] * 0.7),1,n,m, Neurona.Funcion.sigmoidal);
		}
	}

	//endregion

	//region API de la DLL

	/** 
	 Procesar un Vector de Datos en una Red Dada (1: Serie Mensual, 2: Serie Semanal, 3: Serie Mensual)
	 
	 @param numRed net index 
	 @param datosEntrada input data array 
	 @return  mean squared error 
	*/
	public final double procesarVector(int numRed, double[] datosEntrada)
	{
		Red[numRed].procesar(datosEntrada);
		double sal = Red[numRed].getSalida()[0];
		System.Diagnostics.Debug.Write("Entrada (" + datosEntrada.length + ") : ");
		for (int i = 0;i < datosEntrada.length; i++)
		{
			System.Diagnostics.Debug.Write(datosEntrada[i] + ",");
		}
		System.Diagnostics.Debug.WriteLine("\nSalida: " + sal * 100.00);
		return sal * 100.00;

	}

	/** 
	 Procesar un Archivo ASCII con salida en Log
	 
	 @param numRed net index 
	 @param path file path 
	*/
	public final void procesarArchivo(int numRed, String path)
	{

		double umbral = 0.5;
		ArchivoLectura arch = new ArchivoLectura(path,Size[numRed],'\t');
		ArrayList Datos = arch.generarVector(Size[numRed]);

//C# TO JAVA CONVERTER WARNING: The java.io.OutputStreamWriter constructor does not accept all the arguments passed to the System.IO.StreamWriter constructor:
//ORIGINAL LINE: StreamWriter salida = new StreamWriter("..\\..\\..\\out.txt",false);
		java.io.OutputStreamWriter salida = new java.io.OutputStreamWriter("..\\..\\..\\out.txt");
		salida.write("PROCESO DE DATOS : PREDECIBILIDAD DE SERIES NUMERICAS" + System.lineSeparator());
		salida.write("" + System.lineSeparator());
		salida.write("" + System.lineSeparator());

		salida.write("NUM\t" + "% PRED." + "\t" + "JUICIO" + System.lineSeparator());
		salida.write("" + System.lineSeparator());

		int i = 1;
		for (double[] datosEntrada : Datos)
		{
			double sal = this.procesarVector(numRed, datosEntrada);
			salida.write(i++ + ")" + "\t" + (new Double(sal)).toString("//0.//0") + "% " + "\t" + ((sal > umbral)?"1":"0") + System.lineSeparator());
		}
		salida.write("" + System.lineSeparator());
		salida.close();
	}

	/** 
	 Entrenar la Red 
	 
	 @param numRed net index 
	 @param Datos trainning set 
	 @param ciclos number of epochs 
	 @return  mean squared error 
	*/
	public final double entrenar(int numRed, ArrayList Datos, int ciclos)
	{

		double mse = 0.00;
		for (int c = 0;c < ciclos;c++)
		{
			for (double[] datosEntrada : Datos)
			{
				Red[numRed].entrenar(datosEntrada);
			}
		}
		for (double[] datosEntrada : Datos)
		{
			mse += Red[numRed].entrenar(datosEntrada);
		}
		return mse / Datos.size();
	}

	/** 
	 Entrenar la Red en base a un Archivo ASCII previamente seteado un determinado numero de ciclos
	 
	 @param numRed net index 
	 @param path file path 
	 @param ciclos number of epochs 
	 @return  mean squared error 
	*/
	public final double entrenarDeArchivo(int numRed, String path, int ciclos)
	{

		ArchivoLectura arch = new ArchivoLectura(path,Size[numRed] + 1,'\t');
		ArrayList Datos = arch.generarVector(Size[numRed] + 1);

		double mse = 0.00;
		for (int c = 0;c < ciclos;c++)
		{
			for (double[] datosEntrada : Datos)
			{
				Red[numRed].entrenar(datosEntrada);
			}
		}
		for (double[] datosEntrada : Datos)
		{
			mse += Red[numRed].entrenar(datosEntrada);
		}
		return mse / Datos.size();
	}

	/** 
	 Setter del Coeficiente de Entrenamiento N para una determinada red (numRed).
	 
	 @param numRed net number 
	 @param valor trainning coefficent value 
	*/
	public final void setCoefEntrenamiento(int numRed, double valor)
	{
		this.Red[numRed].N = valor;
	}


	/** 
	 Guardar una Configuraci�n de pesos de una determinada red (1: Mensual, 2: Semanal, 3: Diaria) en un stream dado.
	 
	 @param numRed net number 
	 @param pesos stream with weights 
	*/
	public final void saveTo(int numRed, java.io.FileOutputStream pesos)
	{
		Red red = Red[numRed];
		String linea;
		Neurona ne;
		int indice;
		for (int i = 0;i < red.Oculta;i++)
		{
			ne = red.getNeurona(1,i);
			linea = 1 + " " + i + " " + red.Entrada + "";
			for (indice = 0;indice < linea.length();indice++)
			{
//C# TO JAVA CONVERTER WARNING: Unsigned integer types have no direct equivalent in Java:
//ORIGINAL LINE: pesos.WriteByte((byte)linea[indice]);
				pesos.write((byte)linea.charAt(indice));
			}
//C# TO JAVA CONVERTER WARNING: Unsigned integer types have no direct equivalent in Java:
//ORIGINAL LINE: pesos.WriteByte((byte)' ');
			pesos.write((byte)' ');
			for (int j = 0;j < red.Entrada;j++)
			{
				linea = ne.getW(j) + "";
				for (indice = 0;indice < linea.length();indice++)
				{
//C# TO JAVA CONVERTER WARNING: Unsigned integer types have no direct equivalent in Java:
//ORIGINAL LINE: pesos.WriteByte((byte)linea[indice]);
					pesos.write((byte)linea.charAt(indice));
				}
//C# TO JAVA CONVERTER WARNING: Unsigned integer types have no direct equivalent in Java:
//ORIGINAL LINE: pesos.WriteByte((byte)' ');
				pesos.write((byte)' ');
			}
//C# TO JAVA CONVERTER WARNING: Unsigned integer types have no direct equivalent in Java:
//ORIGINAL LINE: pesos.WriteByte((byte)'\r');
			pesos.write((byte)'\r');
//C# TO JAVA CONVERTER WARNING: Unsigned integer types have no direct equivalent in Java:
//ORIGINAL LINE: pesos.WriteByte((byte)'\n');
			pesos.write((byte)'\n');
		}

		for (int i = 0;i < red.Salida;i++)
		{
			ne = red.getNeurona(2,i);
			linea = "2" + " " + i + " " + red.Oculta + " ";
			for (indice = 0;indice < linea.length();indice++)
			{
//C# TO JAVA CONVERTER WARNING: Unsigned integer types have no direct equivalent in Java:
//ORIGINAL LINE: pesos.WriteByte((byte)linea[indice]);
				pesos.write((byte)linea.charAt(indice));
			}
			for (int j = 0;j < red.Oculta;j++)
			{
				linea = ne.getW(j) + " ";
				for (indice = 0;indice < linea.length();indice++)
				{
//C# TO JAVA CONVERTER WARNING: Unsigned integer types have no direct equivalent in Java:
//ORIGINAL LINE: pesos.WriteByte((byte)linea[indice]);
					pesos.write((byte)linea.charAt(indice));
				}
			}
		}
//C# TO JAVA CONVERTER WARNING: Unsigned integer types have no direct equivalent in Java:
//ORIGINAL LINE: pesos.WriteByte((byte)'f');
		pesos.write((byte)'f');
//C# TO JAVA CONVERTER WARNING: Unsigned integer types have no direct equivalent in Java:
//ORIGINAL LINE: pesos.WriteByte((byte)'\r');
		pesos.write((byte)'\r');
//C# TO JAVA CONVERTER WARNING: Unsigned integer types have no direct equivalent in Java:
//ORIGINAL LINE: pesos.WriteByte((byte)'\n');
		pesos.write((byte)'\n');
	}

	/** 
	 Recuperar una Configuraci�n de pesos optimizada para una determinada Red (1: Mensual, 2: Semanal, 3:Diaria) desde un Stream dado.
	 
	 @param numRed net number 
	 @param pesos stream with weights 
	*/
	public final void readFrom(int numRed, java.io.FileInputStream pesos)
	{
		Red red = Red[numRed];
		char caracter;
		ArrayList lineas = new ArrayList();
		String linea = "";
		int byteInt;
//C# TO JAVA CONVERTER WARNING: Unsigned integer types have no direct equivalent in Java:
//ORIGINAL LINE: int fin = Convert.ToByte('f');
		int fin = (byte)'f';
		while (true)
		{
			byteInt = pesos.read();
			if (byteInt == fin)
			{
				//Debug.WriteLine("[" + linea + "]");
				lineas.add(linea);
				break;
			}
			caracter = (char)byteInt;
			if (caracter == '\r')
			{
				pesos.read();
				//Debug.WriteLine("[" + linea + "]");
				lineas.add(linea);
				linea = "";
			}
			else
			{
				linea += caracter;
			}
		}
		pesos.close();

		char[] separador = {' '};
		String[] valores;
		for (String lin : lineas)
		{
			valores = lin.split(java.util.regex.Pattern.quote(separador.toString()), -1);
			int capa = Integer.parseInt(valores[0]);
			int indice = Integer.parseInt(valores[1]);
			int cantConexiones = Integer.parseInt(valores[2]);
			for (int i = 0;i < cantConexiones;i++)
			{
				red.getNeurona(capa,indice).setW(Double.parseDouble(valores[i + 3]),i);
			}
		}
	}

	/** 
	 Guardar una determinada Configuraci�n de pesos en el archivo de texto pesos+numRed.txt
	 
	 @param numRed net number
	*/
	public final void saveTxt(int numRed)
	{
		Red red = this.Red[numRed];
//C# TO JAVA CONVERTER WARNING: The java.io.OutputStreamWriter constructor does not accept all the arguments passed to the System.IO.StreamWriter constructor:
//ORIGINAL LINE: StreamWriter pesos = new StreamWriter("..\\..\\..\\pesos"+numRed+".txt",false);
		java.io.OutputStreamWriter pesos = new java.io.OutputStreamWriter("..\\..\\..\\pesos" + numRed + ".txt");
		Neurona ne;
		for (int i = 0;i < red.Oculta;i++)
		{
			ne = red.getNeurona(1,i);
			pesos.write("1" + " " + i + " " + red.Entrada + " ");
			for (int j = 0;j < red.Entrada;j++)
			{
				pesos.write(ne.getW(j) + " ");
			}
			pesos.write("" + System.lineSeparator());
		}

		for (int i = 0;i < red.Salida;i++)
		{
			ne = red.getNeurona(2,i);
			pesos.write("2" + " " + i + " " + red.Oculta + " ");
			for (int j = 0;j < red.Oculta;j++)
			{
				pesos.write(ne.getW(j) + " ");
			}
			pesos.write("" + System.lineSeparator());
		}
		pesos.close();
	}

	/** 
	 Recuperar una configuraci�n de pesos optimizada para una red desde el archivo de texto pesos+numRed+.txt
	 
	 @param numRed net number 
	*/
	public final void readTxt(int numRed)
	{
		Red red = Red[numRed];
		java.io.InputStreamReader pesos = new java.io.InputStreamReader("..\\..\\..\\pesos" + numRed + ".txt");
		String linea;
		char[] separador = {' '};
		String[] valores;
		linea = pesos.ReadLine();
		while (linea != null)
		{
			valores = linea.split(java.util.regex.Pattern.quote(separador.toString()), -1);
			int capa = Integer.parseInt(valores[0]);
			int indice = Integer.parseInt(valores[1]);
			int cantConexiones = Integer.parseInt(valores[2]);
			for (int i = 0;i < cantConexiones;i++)
			{
				red.getNeurona(capa,indice).setW(Double.parseDouble(valores[i + 3]),i);
			}
			linea = pesos.ReadLine();
		}
		pesos.close();
	}

	//endregion
}