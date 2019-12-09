package Clustering;

import Clustering.Properties.*;
import NumCalc.*;
import ClassicalStat.*;
import java.util.*;

//region Imports


//endregion


/**  Class Dendrograma for hierarchical one dimension clustering 
*/
public class Dendrogram
{

	//region Fields

	private List<Datum> data;
	private HashMap<Integer, Datum> roots;
	private int count;
	private int counter;
	private HypoTest hypoTest;
	private double significance;
	private int minCount;
	private ECCClustering clust;

	//endregion

	//region Constructor

	/**  Constructor 
	*/
	public Dendrogram()
	{
		roots = new HashMap<Integer, Datum>();
		significance = 0.05;
		minCount = 10;
		hypoTest = new HypoTest(significance, minCount);
		clust = new ECCClustering();
	}

	//endregion

	//region Properties

	/**  Test significance 
	*/
	public final double getSignificance()
	{
		return significance;
	}
	public final void setSignificance(double value)
	{
		significance = value;
		hypoTest = new HypoTest(significance, minCount);
	}

	/**  minimum quantity for testing 
	*/
	public final int getMinCount()
	{
		return minCount;
	}
	public final void setMinCount(int value)
	{
		minCount = value;
		hypoTest = new HypoTest(significance, minCount);
	}

	//endregion

	//region Public Methods


	/**  Calculate from a list of tuples 
	 @param data data for clustering 
	 @param cutoff cut-off 
	 @param homog if homogeneity test for clusters should be performed 
	 @param mixed if mixed kmc should be included 
	 @return  list of lists of data (clusters) 
	*/
	public final ArrayList<ArrayList<Datum>> Calculate(List<Datum> data, int cutoff, boolean homog, boolean mixed)
	{
		if (mixed)
		{
			return CalculateMixed(data, cutoff, homog);
		}
		else
		{
			return CalculatePure(data, cutoff, homog);
		}
	}

	/**  Calculate from a list of values 
	 @param values values for clustering 
	 @param cutoff cut-off 
	 @param homog if homogeneity test for clusters should be performed 
	 @return  list of lists of data (clusters) 
	*/
	public final ArrayList<ArrayList<Datum>> Calculate(List<Double> values, int cutoff, boolean homog)
	{
		LoadData(values);
		ArrayList<ArrayList<Datum>> clusters = Calculate(cutoff, homog);
		return clusters;
	}

	//endregion

	//region Private Methods

	//region Load Data

	private void LoadData(List<Double> values)
	{
		ArrayList<Datum> data = new ArrayList<Datum>();
		Datum d;
		for (int i = 0;i < values.size();i++)
		{
			d = new Datum();
			d.IsLeave = true;
			d.index = i;
			d.level = 0;
			d.value = values.get(i);
			d.left = null;
			d.right = null;
			data.add(d);
		}
		LoadData(data);
	}

	private void LoadData(List<Datum> data)
	{
		this.data = data;
		this.count = data.size();
		this.roots.clear();
		for (Datum d : data)
		{
			d.level = data.size();
			this.roots.put(d.index, d);
		}
		counter = data.get(data.size() - 1).index;
	}

	//endregion 

	//region Calculate

	// region Main Methods

	private ArrayList<ArrayList<Datum>> CalculatePure(List<Datum> data, int cutoff, boolean homog)
	{
		LoadData(data);
		ArrayList<ArrayList<Datum>> clusters = Calculate(cutoff, homog);
		return clusters;
	}

	private ArrayList<ArrayList<Datum>> CalculateMixed(List<Datum> dataList, int cutoff, boolean homog)
	{

		//generate k-means data and single dimension
		Dato[] data = new Dato[dataList.size()];
		Dato dato;
		double[] valores;
		double min = Double.MAX_VALUE;
		double max = -Double.MAX_VALUE;
		for (int i = 0;i < dataList.size();i++)
		{
			valores = new double[1];
			valores[0] = dataList.get(i).value;
			dato = new Dato(valores);
			dato.setCodigo(dataList.get(i).code);
			dato.setDim(1);
			data[i] = dato;
			if (dato.getValores()[0] < min)
			{
				min = dato.getValores()[0];
			}
			if (dato.getValores()[0] > max)
			{
				max = dato.getValores()[0];
			}
		}

		Dimension[] dims = new Dimension[1];
		Dimension dim = new Dimension();
		dim.min = min;
		dim.max = max;
		dim.nDivisions = 20;
		dim.normalize = false;
		dim.logarithm = false;
		dims[0] = dim;

		//perform k-means clustering
		clust.LoadData(dims, data);
		clust.Clustering();

		//generate dendrogram data from clusters (index is position in the array of clusters)
		Datum d;
		ArrayList<Datum> datA = new ArrayList<Datum>();
		for (int i = 0;i < clust.Clusters.getLength();i++)
		{
			d = new Datum();
			d.code = "cluster";
			d.index = i;
			d.value = clust.Clusters[i].Centro[0];
			d.IsLeave = true;
			datA.add(d);
		}

		//build dendrogram
		this.LoadData(datA);
		ArrayList<ArrayList<Datum>> clustersOfClusters = Calculate(cutoff, homog);

		//clasify elements in final clusters
		ArrayList<ArrayList<Datum>> clusters = new ArrayList<ArrayList<Datum>>();
		ArrayList<Datum> cluster;
		for (ArrayList<Datum> clOfCl : clustersOfClusters)
		{
			cluster = new ArrayList<Datum>();
			Datum dat;
			for (Datum da : clOfCl)
			{
				for (Dato daa : clust.Clusters[da.index].Datos)
				{
					dat = new Datum();
					dat.code = daa.getCodigo();
					dat.value = daa.getValores()[0];
					cluster.add(dat);
				}
			}
			if (cluster.size() > 0)
			{
				clusters.add(cluster);
			}
		}
		return clusters;
	}

	//endregion

	//region Auxiliar Methods

	private ArrayList<ArrayList<Datum>> Calculate(int cutoff, boolean homog)
	{
		ArrayList<ArrayList<Datum>> clusters = new ArrayList<ArrayList<Datum>>();
		ArrayList<Datum> cluster;
		if (cutoff > roots.size())
		{
			throw new RuntimeException(Strings.Cut_off_must_be_less_or_equal_than_data_count);
		}
		if (!homog && cutoff == roots.size())
		{
			for (Datum d : roots.values())
			{
				cluster = new ArrayList<Datum>();
				cluster.add(d);
				clusters.add(cluster);
			}
			return clusters;
		}

		return BuildDendrogram(roots, cutoff, homog);
	}

	private ArrayList<ArrayList<Datum>> BuildDendrogram(HashMap<Integer, Datum> roots, int cutoff, boolean homog)
	{
		Datum d1 = new Datum();
		Datum d2 = new Datum();
		while (true)
		{
			tangible.RefObject<Clustering.Dendrogram.Datum> tempRef_d1 = new tangible.RefObject<Clustering.Dendrogram.Datum>(d1);
			tangible.RefObject<Clustering.Dendrogram.Datum> tempRef_d2 = new tangible.RefObject<Clustering.Dendrogram.Datum>(d2);
			FindNearest(roots, tempRef_d1, tempRef_d2);
		d2 = tempRef_d2.argValue;
		d1 = tempRef_d1.argValue;
			if (roots.size() == 1 || (homog && roots.size() <= cutoff && !AreHomogeneous(d1, d2)))
			{
				return GetClusters(roots);
			}

			Datum joined = Join(d1, d2);
			joined.level = Math.max(d1.level, d2.level) + 1;
			roots.remove(d1.index);
			roots.remove(d2.index);
			roots.put(joined.index, joined);
			if (!homog && roots.size() == cutoff)
			{
				return GetClusters(roots);
			}
		}
	}

	private void FindNearest(HashMap<Integer, Datum> roots, tangible.RefObject<Datum> datum1, tangible.RefObject<Datum> datum2)
	{
		double dist;
		double minDist = Double.MAX_VALUE;
		for (Datum d1 : roots.values())
		{
			for (Datum d2 : roots.values())
			{
				dist = EuclideanDist(d1.value, d2.value);
				if (d1.index != d2.index && dist < minDist)
				{
					minDist = dist;
					datum1.argValue = d1;
					datum2.argValue = d2;
				}
			}
		}
	}

	private double EuclideanDist(double d1, double d2)
	{
		return Math.sqrt(Math.pow((d1 - d2), 2));
	}

	private Datum Join(Datum d1, Datum d2)
	{
		Datum joined = new Datum();
		joined.index = ++counter;
		joined.value = (d1.value + d2.value) / 2.0;
		joined.left = d1;
		joined.right = d2;
		d1.root = joined;
		d2.root = joined;
		return joined;
	}

	private void GetSubTreeData(Datum d, ArrayList<Datum> data)
	{
		if (d == null)
		{
			return;
		}

		if (d.IsLeave)
		{
			data.add(d);
		}
		else
		{
			GetSubTreeData(d.left, data);
			GetSubTreeData(d.right, data);
		}
	}

	private boolean AreHomogeneous(Datum d1, Datum d2)
	{
		ArrayList<Datum> data1 = new ArrayList<Datum>();
		GetSubTreeData(d1, data1);
		ArrayList<Double> values1 = new ArrayList<Double>();
		for (Datum d : data1)
		{
			values1.add(d.value);
		}

		ArrayList<Datum> data2 = new ArrayList<Datum>();
		GetSubTreeData(d2, data2);
		ArrayList<Double> values2 = new ArrayList<Double>();
		for (Datum d : data2)
		{
			values2.add(d.value);
		}

		double pValue = hypoTest.WaldWolfowitz(values1, values2, true);
		return (pValue > hypoTest.Significance);
	}

	private ArrayList<ArrayList<Datum>> GetClusters(HashMap<Integer, Datum> roots)
	{
		ArrayList<ArrayList<Datum>> clusters = new ArrayList<ArrayList<Datum>>();
		ArrayList<Datum> cluster;
		for (Datum d : roots.values())
		{
			cluster = new ArrayList<Datum>();
			GetSubTreeData(d, cluster);
			clusters.add(cluster);
		}
		return clusters;
	}

	//endregion

	//endregion

	//endregion

	//region Class Datum

	/**  Class Datum (node) 
	*/
	public static class Datum
	{

		/**  Constructor 
		*/
		public Datum()
		{
			left = null;
			right = null;
			IsLeave = false;
		}

		/**  datum index for dictionary 
		*/
		public int index;
		/**  datum code 
		*/
		public String code;
		/**  datum level 
		*/
		public int level;
		/**  datum value 
		*/
		public double value;
		/**  datum root, in the binary tree 
		*/
		public Datum root;
		/**  datum left branch , in the binary tree 
		*/
		public Datum left;
		/**  datum right branch , in the binary tree 
		*/
		public Datum right;
		/**  if the datum is a leave (original datum) 
		*/
		public boolean IsLeave;
	}

	//endregion
}