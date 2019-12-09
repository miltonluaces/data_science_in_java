package Clustering;

import NumCalc.*;
import ClassicalStat.*;
import java.util.*;

//region Imports


//endregion


/**  Algorithms repository for K-Means Clustering 
*/

public class KMC
{

	//region Fields 

	private Dimension[] dim;
	private Dato[] datos;
	private int nDim;
	private int nDatos;
	private Cluster[] clusters;
	private int nClusters;
	private double difCentros;
	private RandomGen randGen;
	private Functions func;

	//endregion

	//region Constructors 

	/**  Constructor 
	*/
	public KMC()
	{
		this.func = new Functions();
		this.randGen = new RandomGen();
	}

	/**  Constructor 
	 @param dim array of dimensions 
	 @param datos array of data 
	*/
	public KMC(Dimension[] dim, Dato[] datos)
	{
		this.nDatos = datos.length;
		this.nDim = dim.length;
		this.dim = dim;
		this.datos = datos;
		nClusters = 1;
		for (int i = 0;i < nDim;i++)
		{
			nClusters *= dim[i].nDivisions;
		}
		clusters = new Cluster[nClusters];
		this.randGen = new RandomGen();
		this.func = new Functions();
	}

	/** 
	 Constructor for category dimensions
	 
	 @param scalDim scalar dimensions
	 @param scalDatos scalar data
	 @param catDim categoric dimensions
	 @param catDatos categoric data
	*/
	public KMC(Dimension[] scalDim, Dato[] scalDatos, CatDimension[] catDim, ArrayList<ArrayList<boolean[]>> catDatos)
	{
		clusters = new Cluster[nClusters];
		this.randGen = new RandomGen();
		this.func = new Functions();

		this.nDim = scalDim.length;
		for (int i = 0; i < catDim.length; i++)
		{
		   this.nDim += catDim[i].nCats;
		}

		this.dim = new Dimension[nDim];
		for (int i = 0; i < scalDim.length; i++)
		{
			this.dim[i] = scalDim[i];
		}
		Dimension cd;
		int ind = scalDim.length;
		for (int i = 0; i < catDim.length; i++)
		{
			for (int j = 0; j < catDim[i].nCats;j++)
			{
				cd = new Dimension();
				cd.index = ind++;
				cd.logarithm = false;
				cd.normalize = false;
				cd.min = 0;
				cd.max = 1;
				cd.nDivisions = 2;
				cd.sort = false;
				cd.weight = catDim[i].weight;
				this.dim[cd.index] = cd;
			}
		}

		this.nDatos = scalDatos.length;
		this.datos = new Dato[nDatos];

		Dato d;
		for (int i = 0; i < nDatos; i++)
		{
			//foreach scalDim a value from the same scalDato to d
			d = new Dato(this.nDim);
			datos[i] = d;
			for (int j = 0; j < scalDim.length; j++)
			{
				d.Valores[j] = scalDatos[i].Valores[j];
				if (d.Valores[j] < dim[j].min)
				{
					dim[j].min = d.Valores[j];
				}
				if (d.Valores[j] > dim[j].max)
				{
					dim[j].max = d.Valores[j];
				}
			}
			//foreach each value of each catDim add to d a value from the correspondant catDato 
			for (int index = scalDim.length; index < nDim; index++)
			{
				for (int j = 0; j < catDim.length; j++)
				{
					for (int k = 0; k < catDim[j].nCats; k++)
					{
						d.Valores[index] = (catDatos.get(i).get(j)[k]) ? 1 : 0;
					}
				}
			}
		}

		nClusters = 1;
		for (int i = 0; i < nDim; i++)
		{
			nClusters *= this.dim[i].nDivisions;
		}
		clusters = new Cluster[nClusters];
		this.randGen = new RandomGen();
		this.func = new Functions();
	}

	//endregion

	//region Properties

	/** 
	 Dimensions
	*/
	public final Dimension[] getDim()
	{
		return dim;
	}

	/** 
	 Data
	*/
	public final Dato[] getDatos()
	{
		return datos;
	}

	/** 
	 Clusters result of the clustering process
	*/
	public final Cluster[] getClusters()
	{
		return clusters;
	}

	//endregion

	//region Public Methods

	/**  Main calculation method 
	*/
	public final void Clustering()
	{
		Reset();
		if (nDim <= 3)
		{
			UbicarCentrosUniforme();
		}
		else
		{
			UbicarCentrosRandom();
		}

		ArrayList newClusters;
		int it = 1;
		while (difCentros > 0)
		{
			 Clasificar(false);
			 if (it == 10)
			 {
				 //eliminacion de clusters vacios
				 newClusters = new ArrayList();
				 for (Cluster cluster : this.clusters)
				 {
					 if (cluster.NDatos != 0)
					 {
						 newClusters.add(cluster);
					 }
				 }
				 this.nClusters = newClusters.size();
				 this.clusters = new Cluster[this.nClusters];
				 for (int j = 0;j < newClusters.size();j++)
				 {
					 this.clusters[j] = (Cluster)newClusters.get(j);
				 }
			 }
			 //Console.Write((it++) + "\n");
			 CalcularCentros();
		}
	}

	//endregion

	//region Private Methods

	//region Place centres

	private void Reset()
	{
		difCentros = Double.MAX_VALUE;
		for (Dato dato : datos)
		{
			dato.Cluster = null;
		}
	}

	private void UbicarCentrosUniforme()
	{
		Dato centro;
		double[] dist = new double[nDim];
		for (int i = 0;i < nDim;i++)
		{
			dist[i] = (dim[i].max - dim[i].min) / dim[i].nDivisions;
		}

		switch (nDim)
		{
			case 1:
				for (int i = 0;i < dim[0].nDivisions;i++)
				{
					centro = new Dato(nDim);
					centro[0] = dim[0].min + i * dist[0];
					Normalize(centro);
					clusters[i] = new Cluster(centro); //Console.WriteLine(centro);
				}
				break;

			case 2:
				int index = 0;
				for (int i = 0;i < dim[0].nDivisions;i++)
				{
					for (int j = 0;j < dim[1].nDivisions;j++)
					{
						centro = new Dato(nDim);
						centro[0] = dim[0].min + i * dist[0];
						centro[1] = dim[1].min + j * dist[1];
						Normalize(centro);
						//Console.WriteLine("c0 = " + centro[0] + ", cn0 = " + centro.Norm[0] + " ; " + " c1 = " + centro[1] + ", cn1 = " + centro.Norm[1]);
						clusters[index++] = new Cluster(centro); //Console.WriteLine(index + ") "+ centro);
					}
				}
				break;

			case 3:
				index = 0;
				for (int i = 0;i < dim[0].nDivisions;i++)
				{
					for (int j = 0;j < dim[1].nDivisions;j++)
					{
						for (int k = 0;k < dim[2].nDivisions;k++)
						{
							centro = new Dato(nDim);
							centro[0] = dim[0].min + i * dist[0];
							centro[1] = dim[1].min + j * dist[1];
							centro[2] = dim[2].min + k * dist[2];
							Normalize(centro);
							//Console.WriteLine("c0 = " + centro[0] + ", cn0 = " + centro.Norm[0] + " ; " + " c1 = " + centro[1] + ", cn1 = " + centro.Norm[1]);
							clusters[index++] = new Cluster(centro); //Console.WriteLine(index + ") "+ centro);
						}
					}
				}
				break;

		}
	}

	private void UbicarCentrosRandom()
	{
		Dato centro;
		randGen.Reset();
		System.Diagnostics.Debug.WriteLine("Centros:");
		for (int i = 0;i < clusters.length;i++)
		{
			centro = new Dato(nDim);
			Normalize(centro);
			for (int j = 0;j < centro.Dim;j++)
			{
				centro[j] = randGen.NextDouble(dim[j].min, dim[j].max);
			}
			clusters[i] = new Cluster(centro);
			System.Diagnostics.Debug.WriteLine(centro);
		}
	}

	private void UbicarCentrosPorDensidad(double radio)
	{

		//Seteo de la densidad
		Dato maxDato = null;
		double dens;
		double maxDensity = 0.0;
		for (Dato dato : datos)
		{
			dato.Densidad = 0.0;
			for (Dato otroDato : datos)
			{
				if (otroDato != dato)
				{
					dens = 1 / func.DistEuclid(dato.Valores, otroDato.Valores);
					dato.Densidad += dens;
					System.Diagnostics.Debug.Write((new Double(dens)).toString("//0.000") + " ");
				}
			}
			dato.Densidad /= datos.length;
			System.Diagnostics.Debug.WriteLine("\n Dato: " + dato);
			if (dato.Densidad >= maxDensity)
			{
				maxDensity = dato.Densidad;
				maxDato = dato;
			}
		}

		//Ubicacion de los Centros
		ArrayList datosDens = new ArrayList();
		datosDens.addAll(Arrays.asList(datos));

		Dato centro;
		System.Diagnostics.Debug.WriteLine("Centros:");
		for (int i = 0;i < clusters.length;i++)
		{
			centro = maxDato;
			clusters[i] = new Cluster(centro);
			System.Diagnostics.Debug.WriteLine(centro);

			maxDensity = 0.0;
			for (Dato dato : datosDens)
			{
				if (func.DistEuclid(centro.Valores,dato.Valores) <= radio)
				{
					datosDens.remove(dato);
				}
				else
				{
					if (dato.Densidad >= maxDensity)
					{
						maxDensity = dato.Densidad;
						maxDato = dato;
					}
				}
			}
		}
	}


	//endregion

	//region Classify and calculate centres

	private void Clasificar(boolean kdt)
	{
		KDTree kdTree = null;
		HashMap<String, Cluster> dicCentroCluster = null;
		double minDist;
		double dist;
		Cluster clusterMasCercano;
		for (int i = 0;i < clusters.length;i++)
		{
			clusters[i].Reset();
		}
		if (kdt)
		{
			dicCentroCluster = new HashMap<String, Cluster>();
			ArrayList<Cluster> clustersList = new ArrayList<Cluster>();
			clustersList.addAll(Arrays.asList(clusters));
			for (Cluster cluster : clusters)
			{
				cluster.Centro.Codigo = cluster.Centro.GetCode();
				dicCentroCluster.put(cluster.Centro.Codigo, cluster);
			}
			kdTree = new KDTree(nDim);
			kdTree.LoadClusters(clustersList);
		}

		for (int i = 0;i < datos.length;i++)
		{
			minDist = Double.MAX_VALUE;
			if (kdt)
			{
				//busqueda por KDTree Nearest Neighbour Search
				Dato centro = kdTree.NearestNeighbourSearch(datos[i]).dato;
				clusterMasCercano = dicCentroCluster.get(centro.Codigo);
			}
			else
			{
				//busqueda por fuerza bruta
				clusterMasCercano = clusters[0];
				for (int j = 0;j < clusters.length;j++)
				{
					dist = func.DistEuclid(datos[i].Valores, clusters[j].Centro.Valores);
					if (dist <= minDist)
					{ // clasifica en el ultimo si equidistantes
						minDist = dist;
						clusterMasCercano = clusters[j];
					}
				}
			}
			clusterMasCercano.Add(datos[i]);
			datos[i].Cluster = clusterMasCercano;
		}
	}

	private void CalcularCentros()
	{
		difCentros = 0.0;
		Dato centro, centroAnt;
		double dif;
		for (int i = 0;i < clusters.length;i++)
		{
			centroAnt = clusters[i].Centro;
			centro = clusters[i].CalcularCentro();
			Normalize(centro);
			dif = func.DistEuclid(centro.Valores, centroAnt.Valores);
			difCentros += dif;
		}
	}

	//endregion

	//region Helper Functions

	private void Normalize(Dato dato)
	{
		for (int i = 0;i < dim.length;i++)
		{
			//no se logaritma ni se normaliza
			if (!dim[i].logarithm && !dim[i].normalize)
			{
				dato.Norm[i] = dato.Valores[i];
			}
			//se logaritma pero no se normaliza
			else if (dim[i].logarithm && !dim[i].normalize)
			{
				dato.Norm[i] = Math.log(dato.Valores[i]);
			}
			//no se logaritma pero se normaliza
			else if (!dim[i].logarithm && dim[i].normalize)
			{
				dato.Norm[i] = (dato.Valores[i] - dim[i].min) / (dim[i].max - dim[i].min);
			}
			//se logaritma y se normaliza
			else if (dim[i].logarithm && dim[i].normalize)
			{
				dato.Norm[i] = (Math.log(dato.Valores[i]) - Math.log(dim[i].min)) / (Math.log(dim[i].max) - (Math.log(dim[i].min)));
			}
		}
	}


	//endregion

	//endregion

	//region Debugging 

	/**  Generate log 
	 @param nClusters number of clusters 
	 @param pond if weights are applied 
	*/
	public final void ClusteringLog(int nClusters, boolean pond)
	{

		//reset
		clusters = new Cluster[nClusters];
		difCentros = Double.MAX_VALUE;
		for (Dato dato : datos)
		{
			dato.Cluster = null;
		}

		//definir ubicacion de los centros
		Dato centro;
		for (int c = 0;c < clusters.length;c++)
		{
			centro = new Dato(nDim);
			Normalize(centro);
			for (int j = 0;j < centro.Dim;j++)
			{
				centro[j] = randGen.NextDouble(dim[j].min, dim[j].max);
				//Console.WriteLine(j+") [" + mins[j] + "," + maxs[j]+ "]");
			}
			clusters[c] = new Cluster(centro);
		}

		//loop principal con impresion
		int i = 1;
		try
		{
//C# TO JAVA CONVERTER WARNING: The java.io.OutputStreamWriter constructor does not accept all the arguments passed to the System.IO.StreamWriter constructor:
//ORIGINAL LINE: StreamWriter so = new StreamWriter("..\\KMC.txt",false);
			java.io.OutputStreamWriter so = new java.io.OutputStreamWriter("..\\KMC.txt");
			so.write("CLUSTERING - K-MEANS : " + datos.length + " Datos con " + dim + " elementos en " + clusters.length + " clusters." + System.lineSeparator());
			so.write("------------------------------------------------------------------------------------------------------" + System.lineSeparator());
			while (difCentros > 0)
			{
				Clasificar(true);
				so.write("\n\nITERACION " + (i++) + " : " + this + System.lineSeparator());
				CalcularCentros();
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

	//region ToString Override 

	/**  ToString override  
	 @return  the formatted string 
	*/
	@Override
	public String toString()
	{
		String ret = String.format(Properties.Strings.Movimientos_centros,difCentros);
		for (int i = 0; i < clusters.length; i++)
		{
			ret = ret + "\n" + Properties.Strings.CLUSTER + (i + 1) + "\n" + clusters[i] + "\n";
		}
		return ret;
	}

	//endregion
}