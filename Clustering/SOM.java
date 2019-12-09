package Clustering;

import NumCalc.*;
import java.util.*;

//region Imports


//endregion


/**  Self-organized map class 
*/
public class SOM
{

	//region Fields

	private int x;
	private int y;
	private Neurona [][] map;
	private int dim;
	private double coefWinner;
	private double coefNeigh;
	private ArrayList datos;
	private double[] max;
	private double[] min;

	//endregion

	//region Constructor

	/**  Constructor 
	 @param dim number of dimensions 
	 @param x projection x axis 
	 @param y projectin y axis 
	 @param coefWinner "winner takes all" coefficent 
	 @param coefNeigh neighbourhood coefficent 
	*/
	public SOM(int dim, int x, int y, double coefWinner, double coefNeigh)
	{
		this.dim = dim;
		this.x = x;
		this.y = y;
		this.coefWinner = coefWinner;
		this.coefNeigh = coefNeigh;
		this.map = new Neurona[x][y];
		for (int i = 0;i < x;i++)
		{
			for (int j = 0;j < y;j++)
			{
				map[i][j] = new Neurona(dim,i,j);
			}
		}
		min = new double[dim];
		max = new double[dim];
		for (int i = 0;i < dim;i++)
		{
			min[i] = Double.MAX_VALUE;
			max[i] = -Double.MAX_VALUE;
		}
	}

	//endregion

	//region Properties, Setters and Getters

	/**  Coordenada x del mapa 
	*/
	public final int getX()
	{
		return x;
	}

	/**  Coordenada y del mapa 
	*/
	public final int getY()
	{
		return y;
	}

	/**  Dimension de los vectores de datos 
	*/
	public final int getDim()
	{
		return dim;
	}

	/**  Coeficiente de entrenamiento ganadora 
	*/
	public final double getCoefWinner()
	{
		return coefWinner;
	}

	/**  Coeficiente de entrenamiento vecindad 
	*/
	public final double getCoefNeigh()
	{
		return coefNeigh;
	}

	/**  Mapa de neuronas 
	*/
	public final Neurona[][] getMap()
	{
		return map;
	}

	/**  Obtener una neurona del mapa 
	 @param x axis x value 
	 @param y axis y value 
	 @return  Neurona of this coordinates 
	*/
	public final Neurona Get(int x, int y)
	{
		return map[x][y];
	}

	/**  Datos para el Clustering 
	*/
	public final ArrayList getDatos()
	{
		return datos;
	}
	public final void setDatos(ArrayList value)
	{
		datos = value;
	}

	//endregion

	//region Public Methods

	//region Load Data

	/**  Load data for clustering  
	 @param datos arraylist of data 
	*/
	public final void LoadData(ArrayList datos)
	{
		this.datos = datos;
		for (Dato dato : datos)
		{
			for (int i = 0;i < dato.Valores.getLength();i++)
			{
				if (dato[i] < min[i])
				{
					min[i] = dato[i];
				}
				if (dato[i] > max[i])
				{
					max[i] = dato[i];
				}
			}
		}
	}

	//endregion

	//region Process 

	/**  Proceso del SOM con Array de Valores 
	 @param valores values for clustering 
	 @return  Active neuron 
	*/
	public final Neurona Process(double[] valores)
	{
		Norma norma = new Norma(valores);
		norma.Normalize(Norma.NormType.minMax);
		double[] valoresNorm = norma.Norm;
		double minError = Double.MAX_VALUE;
		double error = 0.0;
		Neurona winner = null;

		for (int i = 0;i < x;i++)
		{
			for (int j = 0;j < y;j++)
			{
				error = map[i][j].Distancia(valoresNorm); //Console.WriteLine("Dist["+i+"] = "+dist);
				if (error < minError)
				{
					minError = error;
					winner = map[i][j];
				}
			}
		}
		winner.Error = error;
		return winner;
	}

	/**  Proceso del SOM con Datos 
	 @param dato datum for clustering 
	 @return  Active neuron selected 
	*/
	public final Neurona Process(Dato dato)
	{
		return Process(dato.Valores);
	}

	/**  Proceso del SOM con n datos a clasificar 
	 @param datos data for clustering 
	*/
	public final void Process(ArrayList datos)
	{
		this.datos = datos;
		Reset();
		Neurona winner;
		for (Dato dato : datos)
		{
			winner = Process(dato);
			winner.AddDato(dato);
		}
	}

	/**  Procesar datos anteriormente cargados en el SOM 
	*/
	public final void Process()
	{
		Reset();
		Neurona winner;
		for (Dato dato : this.datos)
		{
			winner = Process(dato);
			winner.AddDato(dato);
		}
	}

	//endregion

	//region Entrenamiento

	/**  Entrenamiento del SOM con Array de Valores 
	 @param valores trainning set  
	 @return  mean squared error 
	*/
	public final double Train(double[] valores)
	{
		Norma norma = new Norma(valores);
		norma.Normalize(Norma.NormType.minMax);
		double[] valoresNorm = norma.Norm;
		Neurona winner = Process(valores);
		double totDif = 0.0;
		double dif;

		//entrenamiento de la ganadora
		for (int i = 0;i < dim;i++)
		{
			dif = valoresNorm[i] - winner.GetW(i);
			winner.IncW(i, coefWinner * dif);
			totDif += Math.abs(coefWinner * dif);
		}

		//entrenamiento vecindad
		for (Neurona neigh : GetNeighbors(winner))
		{
			for (int i = 0;i < dim;i++)
			{
				dif = valoresNorm[i] - neigh.GetW(i);
				neigh.IncW(i, coefNeigh * dif);
			}
		}

		return totDif;
	}

	/**  Entrenar con datos anteriormente guardados 
	 @param nCiclos number of epochs 
	*/
	public final void Train(int nCiclos)
	{
		for (int i = 0;i < nCiclos;i++)
		{
			for (Dato dato : this.datos)
			{
				Train(dato.Valores);
			}
		}
	}


	/**  Entrenamiento del SOM n ciclos 
	 @param datos trainning set 
	 @param ciclos number of epochs 
	 @return  mean squared error 
	*/
	public final double Train(double[][] datos, int ciclos)
	{
		double totDif = 0f;
		double dif;
		for (int c = 0;c < ciclos;c++)
		{
			totDif = 0f;
			int cantDatos = datos.length;
			for (int i = 0; i < cantDatos; i++)
			{
				dif = Train(datos[i]);
				totDif += Math.abs(dif);
			}
			//System.out.println("Iteracion "+c+" : diferencia = "+totDif);
		}
		return totDif;
	}

	/**  Automatic trainning method 
	 @param datos data for training (n patterns with n data each) 
	 @return  Mean squared error 
	*/
	public final double AutomaticTrain(double[][] datos)
	{
		double totDif = 0f;
		double dif;
		double totDifAnt = 0f;
		while (totDifAnt != totDif)
		{
			totDifAnt = totDif;
			totDif = 0f;
			int cantDatos = datos.length;
			for (int i = 0; i < cantDatos; i++)
			{
				dif = Train(datos[i]);
				totDif += Math.abs(dif);
			}
			//System.out.println("Iteracion "+(c++)+" : diferencia = "+totDif);
		}
		return totDif;
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
//ORIGINAL LINE: StreamWriter so = new StreamWriter("..\\SOM.txt",false);
			java.io.OutputStreamWriter so = new java.io.OutputStreamWriter("..\\SOM.txt");
			so.write("CLUSTERING - SOM : Datos con " + dim + " elementos en " + (x * y) + " clusters." + System.lineSeparator());
			so.write("------------------------------------------------------------------------------" + System.lineSeparator());
			for (int i = 0;i < x;i++)
			{
				for (int j = 0;j < y;j++)
				{
					so.write(String.valueOf(map[i][j]) + System.lineSeparator());
				}
			}
			so.close();
			System.Diagnostics.Debug.WriteLine("Completado archivo ");
		}
		catch (java.lang.Exception e)
		{
			System.Diagnostics.Debug.WriteLine("Error");
		}
	}

	//endregion

	//endregion

	//region Private Methods

	private ArrayList GetNeighbors(Neurona winner)
	{
		ArrayList neighbors = new ArrayList();
		if (winner.X > 0)
		{
			neighbors.add(map[winner.X - 1][winner.Y]);
		}
		if (winner.X < x - 1)
		{
			neighbors.add(map[winner.X + 1][winner.Y]);
		}
		if (winner.Y > 0)
		{
			neighbors.add(map[winner.X][winner.Y - 1]);
		}
		if (winner.Y < y - 1)
		{
			neighbors.add(map[winner.X][winner.Y + 1]);
		}
		return neighbors;
	}

	private void Reset()
	{
		for (int i = 0;i < x;i++)
		{
			for (int j = 0;j < y;j++)
			{
				map[i][j].Reset();
			}
		}
	}

	//endregion

}