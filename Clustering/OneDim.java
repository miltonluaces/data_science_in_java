package Clustering;

import java.util.*;

/**  Algorithmic repository class for one-dimension clustering based in radium parameters 
*/
public class OneDim
{

	//region Fields

	private ArrayList<Double> datos;
	private ArrayList<Cluster> clusters;
	private double radioInc;
	private double radioClust;
	private Hashtable datoCluster;
	private Hashtable datoIndex;

	//endregion

	//region Constructors

	/**  Constructor 
	 @param datos array of data 
	 @param radioInc increment on radium 
	 @param radioClust cluster radium 
	*/
	public OneDim(ArrayList<Double> datos, double radioInc, double radioClust)
	{
		this.datos = datos;
		this.clusters = new ArrayList<Cluster>();
		this.radioInc = radioInc;
		this.radioClust = radioClust;
		datoCluster = new Hashtable();
		datoIndex = new Hashtable();
		for (int i = 0;i < datos.size();i++)
		{
			datoIndex.put(datos.get(i), i);
		}
	}

	/**  Constructor 
	 @param datos array of data 
	 @param radioInc increment on radium 
	 @param radioClust cluster radium 
	*/
	public OneDim(double[] datos, double radioInc, double radioClust)
	{
		this.datos = new ArrayList<Double>();
		datoCluster = new Hashtable();
		datoIndex = new Hashtable();
		for (int i = 0;i < datos.length;i++)
		{
			this.datos.add(datos[i]);
			datoIndex.put(datos[i], i);
		}
		this.clusters = new ArrayList<Cluster>();
		this.radioInc = radioInc;
		this.radioClust = radioClust;
	}

	//endregion

	//region Properties

	/**  list of clusters 
	*/
	public final ArrayList<Cluster> getClusters()
	{
		return clusters;
	}

	/**  increment in radium 
	*/
	public final double getRadioInc()
	{
		return radioInc;
	}
	public final void setRadioInc(double value)
	{
		radioInc = value;
	}

	/**  cluster radium 
	*/
	public final double getRadioClust()
	{
		return radioClust;
	}
	public final void setRadioClust(double value)
	{
		radioClust = value;
	}



  //endregion

	//region Public Methods

	/**  Main clustering method 
	 @param sort if is must be sorted 
	 @param incBoundaries if boundaries must be included 
	*/
	public final void Clustering(boolean sort, boolean incBoundaries)
	{
		if (sort)
		{
			Collections.sort(datos);
		}
		if (incBoundaries)
		{
			int c = 0;
			double min = datos.get(0);
			Cluster cluster = new Cluster(c++);
			cluster.Add(0, datos.get(0));
			clusters.add(cluster);
			cluster = new Cluster(c++);
			cluster.Add(1, datos.get(1));
			for (int i = 2;i < datos.size() - 1;i++)
			{
				if (Math.abs(datos.get(i) - datos.get(i - 1)) >= radioInc || Math.abs(datos.get(i) - min) > radioClust)
				{
					clusters.add(cluster);
					cluster = new Cluster(c++);
					min = datos.get(i);
				}
				cluster.Add(i, datos.get(i));
				datoCluster.put(datos.get(i), cluster);
			}
			clusters.add(cluster);
			cluster = new Cluster(c++);
			cluster.Add(datos.size() - 1, datos.get(datos.size() - 1));
			clusters.add(cluster);
		}
		else
		{
			int c = 0;
			double min = datos.get(0);
			Cluster cluster = new Cluster(c);
			cluster.Add(0, datos.get(0));
			for (int i = 1;i < datos.size();i++)
			{
				if (Math.abs(datos.get(i) - datos.get(i - 1)) > radioInc || Math.abs(datos.get(i) - min) > radioClust)
				{
					clusters.add(cluster);
					cluster = new Cluster(c++);
					min = datos.get(i);
				}
				cluster.Add(i, datos.get(i));
				datoCluster.put(datos.get(i), cluster);
			}
			clusters.add(cluster);
		}
	}

	/**  Add a new cluster 
	 @param cluster the cluster to add 
	*/
	public final void Add(Cluster cluster)
	{
		clusters.add(cluster);
	}

	/**  Get an array of centres 
	 @return  the array of centres 
	*/
	public final double[] GetCentros()
	{
		double[] centros = new double[clusters.size()];
		for (int i = 0;i < clusters.size();i++)
		{
			clusters.get(i).Calculate();
			centros[i] = clusters.get(i).getCenter();
		}
		return centros;
	}

	/**  Get the cluster to whom a datum belongs 
	 @param dato the datum 
	 @return  its cluster 
	*/
	public final Cluster GetCluster(double dato)
	{
		return (Cluster)datoCluster.get(dato);
	}

	/**  Get a cluster proxy to a certain datum 
	 @param dato the datum 
	 @return  the proxy cluster 
	*/
	public final Cluster GetProxCluster(double dato)
	{
		double minDist = Double.MAX_VALUE;
		Cluster minCluster = null;
		for (Cluster cluster : clusters)
		{
			cluster.Calculate();
			double dist = Math.abs(cluster.getCenter() - dato);
			if (dist < minDist)
			{
				minDist = dist;
				minCluster = cluster;
			}
		}
		return minCluster;
	}

	/**  Get the index of a datum 
	 @param dato the datum 
	 @return  its index 
	*/
	public final int GetIndex(double dato)
	{
		return (int)datoIndex.get(dato);
	}

	//endregion

	//region Inner classes 

	/**  Speciphical Cluster class for one dimension clustering 
	*/
	public static class Cluster
	{

		//region Fields

		private int index;
		private double center;
		private double min;
		private double max;
		private ArrayList<Double> datos;
		private ArrayList<Integer> indexes;


		//endregion

		//region Constructor

		/**  Constructor 
		 @param index index of the cluster 
		*/
		public Cluster(int index)
		{
			this.index = index;
			datos = new ArrayList<Double>();
			indexes = new ArrayList<Integer>();
		}

		//endregion

		//region Properties , Setters and Getters

		/**  Add a new datum 
		 @param index index of the datum 
		 @param dato the datum 
		*/
		public final void Add(int index, double dato)
		{
			datos.add(dato);
			indexes.add(index);
		}

		/**  Index of the cluster 
		*/
		public final int getIndex()
		{
			return index;
		}

		/**  Center of the cluster 
		*/
		public final double getCenter()
		{
			return center;
		}

		/**  Min value of the cluster 
		*/
		public final double getMin()
		{
			return min;
		}

		/**  Max value of the cluster 
		*/
		public final double getMax()
		{
			return max;
		}

		/**  List of data that belong to the cluster 
		*/
		public final ArrayList<Double> getDatos()
		{
			return datos;
		}

		/**  List of indexes of data that belong to the cluster 
		*/
		public final ArrayList<Integer> getIndexes()
		{
			return indexes;
		}

		//endregion  

		//region Public Methods

		/**  Main calculation method 
		*/
		public final void Calculate()
		{
			double tot = 0.0;
			min = Double.MAX_VALUE;
			max = -Double.MAX_VALUE;
			for (double dato : datos)
			{
				tot += dato;
				if (dato < min)
				{
					min = dato;
				}
				if (dato > max)
				{
					max = dato;
				}
			}
			center = tot / datos.size();
		}
		//endregion
	}

	//endregion
}