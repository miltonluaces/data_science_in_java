package Clustering;

import NumCalc.*;
import java.util.*;

//region Imports


//endregion


/**  Algorithms repository for Hierarchical clustering 
*/

public class Hierarchical
{

	//region Fields

	private Dimension[] dim;
	private Dato[] datos;
	private int nDim;
	private int nDatos;
	private ArrayList<Cluster> clusters;
	private Functions func;

	//endregion

	//region Constructor

	/**  Constructor 
	 @param dim array of dimensions 
	 @param datos array of data 
	*/
	public Hierarchical(Dimension[] dim, Dato[] datos)
	{
		this.func = new Functions();
		this.nDatos = datos.length;
		this.nDim = dim.length;
		this.dim = dim;
		this.datos = datos;
		clusters = new ArrayList<Cluster>();
		//inicializacion: todos los datos son clusters
		Cluster cluster;
		for (Dato dato : datos)
		{
			cluster = new Cluster(dato);
			cluster.Add(dato);
			clusters.add(cluster);
		}
	}

	//endregion

	//region Properties

	/**  Clusters, result of the clustering process 
	*/
	public final ArrayList<Cluster> getClusters()
	{
		return clusters;
	}

	//endregion

	//region Public Methods

	/**  Main method 
	 @param nClustersFinal number of final clusters 
	*/
	public final void Clustering(int nClustersFinal)
	{
		if (nClustersFinal >= datos.length)
		{
			return;
		}

		while (clusters.size() > nClustersFinal)
		{
			double[][] distMatrix = GetDistMatrix(clusters);
			Dato coords = GetMinDistCoords(distMatrix);
			Cluster cluster1 = clusters.get((int)coords[0]);
			Cluster cluster2 = clusters.get((int)coords[1]);
			Cluster join = JoinClusters(cluster1, cluster2);
			clusters.add(join);
			clusters.remove(cluster1);
			clusters.remove(cluster2);
		}
	}

	//endregion

	//region Private Methods

	private double[][] GetDistMatrix(ArrayList<Cluster> clusters)
	{
		double[][] distMatrix = new double[clusters.size()][clusters.size()];
		for (int i = 0;i < clusters.size();i++)
		{
			for (int j = 0;j < i;j++)
			{
				distMatrix[i][j] = func.DistEuclid(clusters.get(i).Centro.Valores, clusters.get(j).Centro.Valores);
				distMatrix[j][i] = distMatrix[i][j];
			}
		}
		return distMatrix;
	}

	private Dato GetMinDistCoords(double[][] distMatrix)
	{
		double minDist = Double.MAX_VALUE;
		Dato coords = new Dato(2);
		for (int i = 0;i < distMatrix.length;i++)
		{
			for (int j = 0;j < distMatrix[0].length;j++)
			{
				if (i != j && distMatrix[i][j] < minDist)
				{
					minDist = distMatrix[i][j];
					coords[0] = i;
					coords[1] = j;
				}
			}
		}
		return coords;
	}

	private Cluster JoinClusters(Cluster cluster1, Cluster cluster2)
	{
		Cluster join = new Cluster(cluster1.Centro);
		join.Add(cluster1.Datos);
		join.Add(cluster2.Datos);
		join.CalcularCentro();
		return join;
	}

	//endregion

}