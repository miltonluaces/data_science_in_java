package Clustering;

import NumCalc.*;
import java.util.*;

//region Imports


//endregion


/**  KD-Tree Lazy clustering for updates (classification in sub-spaces) 
*/
public class KDTree
{

	//region Fields

	private Functions func;
	private List<Cluster> clusters;
	private HashMap<String, Cluster> clustersDic;
	private ArrayList<Dato>[] sortedLists;
	private ArrayList<Dato> centros;
	private Nodo raiz;
	private int nDim;

	//endregion

	//region Constructor

	/**  Constructor 
	 @param nDim number of dimensions 
	*/
	public KDTree(int nDim)
	{
		this.func = new Functions();
		this.nDim = nDim;
		this.clustersDic = new HashMap<String, Cluster>();
	}

	//endregion

	//region Properties

	/**  Root of the hierarchy 
	*/
	public final Nodo getRaiz()
	{
		return raiz;
	}

	/**  List of clusters 
	*/
	public final List<Cluster> getClusters()
	{
		return clusters;
	}

	//endregion


	//region Public Methods

	/**  Load data 
	 @param clusters collection of clusters to load 
	*/
	public final void LoadClusters(List<Cluster> clusters)
	{
		this.clusters = clusters;
		centros = new ArrayList<Dato>();
		clustersDic.clear();
		for (Cluster cluster : clusters)
		{
			centros.add(cluster.Centro);
			clustersDic.put(cluster.Centro.GetCode(), cluster);
		}
		sortedLists = new ArrayList<Dato>[nDim];
		for (int i = 0;i < nDim;i++)
		{
			sortedLists[i] = Sort(centros, i);
		}
		Dato mins = new Dato(nDim);
		Dato maxs = new Dato(nDim);
		for (int i = 0;i < nDim;i++)
		{
			mins[i] = 0.0;
			maxs[i] = 100.0;
		}
		raiz = CreateKDT(centros, 0, centros.size() - 1, 0, -1, mins, maxs);
	}

	/** 
	 Main method
	 
	 @param datos list of data 
	*/
	public final void Discriminate(ArrayList<Dato> datos)
	{

		Dato centro;
		for (Dato dato : datos)
		{
			centro = NearestNeighbourSearch(dato).dato;
			GetCluster(centro).Add(dato);
		}
	}

	/**  Main method 
	 @param datos list of data 
	 @param maxStackSize max size of the stack 
	*/
	public final void Discriminate(ArrayList<Dato> datos, int maxStackSize)
	{
		ClustCalculator cc;
		ArrayList<Thread> threads = new ArrayList<Thread>();
		Thread thread;
		for (Dato dato : datos)
		{
			cc = new ClustCalculator(this, dato);
			thread = new Thread(new ThreadStart(cc.Calculate), maxStackSize);
			threads.add(thread);
			thread.start();
		}

		for (Thread th : threads)
		{
			th.join();
		}
	}

	/**  Get the nearest neighbour node to the datum 
	 @param dato the datum 
	 @return  the nearest neighbour 
	*/
	public final Nodo NearestNeighbourSearch(Dato dato)
	{

		//busqueda inicial en el kd-tree
		Nodo centro = KDTSearch(dato);
		if (centro == null)
		{
			return null;
		}

		//chequeo de distancias en la hiperesfera
		double dist = func.DistEuclid(dato.Valores, centro.dato.Valores);
		boolean outOfBounds = false;
		for (int i = 0;i < nDim;i++)
		{
			if (Math.abs(dato.Valores[i] - dist) < centro.mins[i] || Math.abs(dato.Valores[i] + dist) > centro.maxs[i])
			{
				outOfBounds = true;
			}
		}
		if (outOfBounds)
		{
			Dato best = GetBestInSphere(dato, dist);
			if (best.Dim > 0)
			{
				centro = new Nodo(best);
			}
		}
		return centro;
	}

	/**  Get a cluster from its centre 
	 @param centro the centre 
	 @return  the cluster 
	*/
	public final Cluster GetCluster(Dato centro)
	{
		return clustersDic.get(centro.GetCode());
	}

	/**  Clear of data all clusters 
	*/
	public final void ResetClusters()
	{
		for (Cluster cluster : clusters)
		{
			cluster.Reset();
		}
	}

	//endregion

	//region Private Methods

	//region KDTree 

	/**  Search in sub spaces by KD-tree 
	 @param dato the value to search 
	 @return  the nearest neighbour node 
	*/
	public final Nodo KDTSearch(Dato dato)
	{
		Nodo nodo = new Nodo(dato);
		double dist = func.DistEuclid(dato.Valores, raiz.dato.Valores);
		return KDTSearch(nodo, raiz);
	}

	/**  KD-Tree search in sub-spaces, recursive method  
	 @param nodo the node 
	 @param centro actual centre 
	 @return  the node 
	*/
	public final Nodo KDTSearch(Nodo nodo, Nodo centro)
	{
		double dist = func.DistEuclid(nodo.dato.Valores, centro.dato.Valores);
		int d = centro.dim;
		if (nodo.dato.Valores[d] == centro.dato.Valores[d])
		{
			if (nodo.dato == centro.dato)
			{
				return centro;
			}
			else
			{
				if (centro.der == null)
				{
					return centro;
				}
				return KDTSearch(nodo, centro.der);
			}
		}
		else if (nodo.dato.Valores[d] < centro.dato.Valores[d])
		{
			if (centro.izq == null)
			{
				return centro;
			}
			else
			{
				return KDTSearch(nodo, centro.izq);
			}
		}
		else if (nodo.dato.Valores[d] > centro.dato.Valores[d])
		{
			if (centro.der == null)
			{
				return centro;
			}
			else
			{
				return KDTSearch(nodo, centro.der);
			}
		}
		return null;
	}

	/**  Factory method for the KDT 
	 @param datos List of data 
	 @param ini initial 
	 @param fin final 
	 @param h height 
	 @param dim number of dimensions 
	 @param mins datum that contains array of minima dimensions 
	 @param maxs datum that contains array of maxima dimensions 
	 @return  root node 
	*/
	private Nodo CreateKDT(ArrayList<Dato> datos, int ini, int fin, int h, int dim, Dato mins, Dato maxs)
	{
		dim = (dim + 1) % nDim;
		ArrayList<Dato> subDatos = Split(ini, fin, dim);
		Nodo nodo;
		switch (subDatos.size())
		{
			case 0:
				return null;
			case 1:
				nodo = new Nodo(subDatos.get(0), 0, h, dim, mins, maxs);
				nodo.izq = null;
				nodo.der = null;
				return nodo;
			case 2:
				nodo = new Nodo(subDatos.get(0), 0, h, dim, mins, maxs);
				nodo.izq = null;
				Dato minsRama = mins.Clone();
				minsRama[dim] = nodo.dato[dim];
				Nodo rama = new Nodo(subDatos.get(1), 0, h, (dim + 1) % nDim, minsRama, maxs);
				rama.raiz = nodo;
				nodo.der = rama;
				return nodo;
			default:
				int medianIndex = (int)(subDatos.size() / 2);
				Dato median = subDatos.get(medianIndex);
				nodo = new Nodo(median, medianIndex, h, dim, mins, maxs);
				Dato maxsIzq = maxs.Clone();
				maxsIzq[dim] = median[dim];
				Dato minsDer = mins.Clone();
				minsDer[dim] = median[dim];
				nodo.izq = CreateKDT(subDatos, 0, medianIndex - 1, h + 1, dim, mins, maxsIzq);
				nodo.der = CreateKDT(subDatos, medianIndex + 1, subDatos.size() - 1, h + 1, dim, minsDer, maxs);
				return nodo;
		}
	}

	//endregion

	//region Metodos Auxiliares

	private Dato GetBestInSphere(Dato dato, double radio)
	{
		Dato best = new Dato(0);
		double bestDist = radio;
		ArrayList<Dato> otros = new ArrayList<Dato>();
		double dist;
		for (Dato centro : centros)
		{
			dist = func.DistEuclid(centro.Valores, dato.Valores);
			if (dist < bestDist)
			{
				bestDist = dist;
				best = centro;
			}
		}
		return best;
	}

	private ArrayList<Dato> Sort(ArrayList<Dato> datos, int sortDim)
	{

		ArrayList<Dato> sorted = new ArrayList<Dato>();
		Dato clon;
		for (Dato dato : datos)
		{
			clon = dato.Clone();
			clon.SortDim = sortDim;
			sorted.add(clon);
		}
		Collections.sort(sorted);
		return sorted;
	}

	private ArrayList<Dato> Split(int inicio, int fin, int sortDim)
	{

		ArrayList<Dato> splitted = new ArrayList<Dato>();
		for (int i = inicio;i <= fin;i++)
		{
			splitted.add(sortedLists[sortDim].get(i));
		}
		return splitted;
	}

	//endregion

	//endregion

	//region Nodo Class

	/**  Nodo Container Class 
	*/
	public static class Nodo
	{

		/**  the dato itself 
		*/
		public Dato dato;

		/**  wrapper that contains mimima 
		*/
		public Dato mins;

		/**  wrapper that contains maxima 
		*/
		public Dato maxs;

		/**  location parameter 
		*/
		public int loc;

		/**  height parameter 
		*/
		public int h;

		/**  number of dims 
		*/
		public int dim;

		/**  left branch in Kd-tree 
		*/
		public Nodo izq;

		/**  right branch in Kd-tree  
		*/
		public Nodo der;

		/**  root of Kd-tree 
		*/
		public Nodo raiz;

		/**  Constructor 
		 @param dato a datum to construct a unitarian node 
		*/
		public Nodo(Dato dato)
		{
			this.dato = dato;
			this.loc = 0;
			this.h = 0;
			this.dim = 0;
			mins = new Dato(dato.Dim);
			maxs = new Dato(dato.Dim);
			for (int i = 0;i < dato.Dim;i++)
			{
				mins[i] = 0.0;
				maxs[i] = 100.0;
			}
			Nodo nodo = new Nodo(dato, 0, 0, 0, mins, maxs);
		}

		/**  Constructor
		 @param dato a datum to construct a unitarian node 
		 @param loc location 
		 @param h height 
		 @param dim number of dimensions 
		 @param mins datum that contains array of minima dimensions 
		 @param maxs  datum that contains array of maxima dimensions 
		*/
		public Nodo(Dato dato, int loc, int h, int dim, Dato mins, Dato maxs)
		{
			this.dato = dato;
			this.loc = loc;
			this.h = h;
			this.dim = dim;
			this.mins = mins;
			this.maxs = maxs;
		}

		/**  Get value of current dim 
		 @return  value of current dim 
		*/
		public final double GetDimValue()
		{
			return dato.Valores[dim];
		}

		/**  ToString override  
		 @return  the formatted string 
		*/
		@Override
		public String toString()
		{
			return dato.toString();
		}
	}

	//endregion

}