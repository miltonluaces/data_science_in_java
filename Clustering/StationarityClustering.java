package Clustering;

import Clustering.Properties.*;
import ClassicalStat.*;
import java.util.*;

/**  Algorithmic repository class for stationarity time series clustering 
*/
public class StationarityClustering
{

	//region Fields

	private double origMaxDiff;
	private double origMaxStDev;
	private double maxDiff;
	private double maxStDev;
	private double maxStDevDiff;
	private int minCount;

	//endregion

	//region Constructor

	/** 
	 Constructor 
	 
	 @param maxDiff maximum difference allowed 
	 @param maxStDev maximum standard deviation allowed 
	 @param minCount minimum number of elements allowed 
	*/
	public StationarityClustering(double maxDiff, double maxStDev, int minCount)
	{
		this.maxDiff = maxDiff;
		this.maxStDev = maxStDev;
		this.origMaxDiff = maxDiff;
		this.origMaxStDev = maxStDev;
		this.minCount = minCount;
	}

	//endregion

	//region Properties

	/**  maximum difference allowed within a cluster 
	*/
	public final double getMaxDiff()
	{
		return maxDiff;
	}
	public final void setMaxDiff(double value)
	{
		maxDiff = value;
		origMaxDiff = maxDiff;
	}

	/**  maximum standard deviation difference allowed within a cluster 
	*/
	public final double getMaxStDevDiff()
	{
		return maxStDevDiff;
	}
	public final void setMaxStDevDiff(double value)
	{
		maxStDevDiff = value;
	}

	/**  original maximum difference allowed within a cluster 
	*/
	public final double getOrigMaxDiff()
	{
		return origMaxDiff;
	}
	public final void setOrigMaxDiff(double value)
	{
		origMaxDiff = value;
	}

	/**  maximum standard deviation allowed within a cluster 
	*/
	public final double getMaxStDev()
	{
		return maxStDev;
	}
	public final void setMaxStDev(double value)
	{
		maxStDev = value;
		origMaxStDev = maxStDev;
	}

	/**  original maximum standard deviation allowed within a cluster 
	*/
	public final double getOrigMaxStDev()
	{
		return origMaxStDev;
	}
	public final void setOrigMaxStDev(double value)
	{
		origMaxStDev = value;
	}

	/**  minimum number of values allowed within a cluster 
	*/
	public final int getMinCount()
	{
		return minCount;
	}
	public final void setMinCount(int value)
	{
		minCount = value;
	}

	//endregion

	//region Setters and Getters

	/**  Calculate mean from a list of clusters 
	 @param clusters the list of clusters 
	 @return  the calculated mean 
	*/
	public final ArrayList<Double> GetMeans(ArrayList<StatCluster> clusters)
	{
		ArrayList<Double> means = new ArrayList<Double>();
		for (StatCluster sc : clusters)
		{
			for (int i = sc.getFirstIndex();i <= sc.getLastIndex();i++)
			{
				means.add(sc.getMean());
			}
		}
		return means;
	}

	/**  Calculate standard deviation from a list of clusters 
	 @param clusters the list of clusters 
	 @return  the calculated standard deviation 
	*/
	public final ArrayList<Double> GetStDevs(ArrayList<StatCluster> clusters)
	{
		ArrayList<Double> stDevs = new ArrayList<Double>();
		for (StatCluster sc : clusters)
		{
			for (int i = sc.getFirstIndex();i <= sc.getLastIndex();i++)
			{
				stDevs.add(sc.getStDev());
			}
		}
		return stDevs;
	}

	/**  Limit indexes 
	 @param clusters the list of clusters 
	 @return  limit 
	*/
	public final ArrayList<Integer> GetLimitIndexes(ArrayList<StatCluster> clusters)
	{
		ArrayList<Integer> limitIndexes = new ArrayList<Integer>();
		for (StatCluster sc : clusters)
		{
			limitIndexes.add(sc.getFirstIndex());
		}
		return limitIndexes;
	}

	//endregion

	//region Public Methods

	/**  Main calculation method 
	 @param serie list containing the serie 
	 @return  list of clusters 
	*/
	public final ArrayList<StatCluster> Clustering(ArrayList<Double> serie)
	{
		if (serie == null)
		{
			throw new RuntimeException(Strings.Serie_is_null);
		}

		ArrayList<StatCluster> clusters = CreateClusters(serie);
		maxDiff = origMaxDiff;
		maxStDev = origMaxStDev;
		while (maxDiff < 100.0)
		{
			JoinNeighbourClusters(clusters);
			JoinInnerClusters(clusters);
			CloseBigClusters(clusters);
			maxDiff += origMaxDiff;
			maxStDev += origMaxStDev;
		}
		return clusters;
	}

	/**  Method for performing recursive clustering upon clusters 
	 @param clusters list of clusters 
	 @return  resulting clusters 
	*/
	public final ArrayList<StatCluster> RecursiveClustering(ArrayList<StatCluster> clusters)
	{
		ArrayList<StatCluster> recursiveClusters = new ArrayList<StatCluster>();
		recursiveClusters.add(clusters.get(clusters.size() - 1));
		for (int i = 0;i < clusters.size() - 2;i++)
		{
			recursiveClusters.add(new StatCluster(clusters.get(clusters.size() - 2 - i), recursiveClusters.get(i)));
		}
		return recursiveClusters;
	}

	//endregion 

	//region Private Methods

	//region Calculate Methods

	private ArrayList<StatCluster> CreateClusters(ArrayList<Double> serie)
	{
		ArrayList<StatCluster> clusters = new ArrayList<StatCluster>();
		StatCluster sC;
		for (int j = 0;j < serie.size();j++)
		{
			sC = new StatCluster(j, serie.get(j));
			clusters.add(sC);
		}
		return clusters;
	}

	private void JoinNeighbourClusters(ArrayList<StatCluster> clusters)
	{
		int i;
		ArrayList<Integer> joined = new ArrayList<Integer>();
		double diff1, diff2, stDev1, stDev2, stDevDiff1, stDevDiff2, countProp1, countProp2;
		while (true)
		{
			i = 1;
			joined.clear();
			while (i < clusters.size() - 1)
			{
				if (clusters.get(i - 1).getClose() || clusters.get(i).getClose() || clusters.get(i + 1).getClose())
				{
					i = i + 1;
				}
				else
				{
					diff1 = GetDiff(clusters.get(i - 1), clusters.get(i));
					diff2 = GetDiff(clusters.get(i), clusters.get(i + 1));
					stDev1 = GetStDev(clusters.get(i - 1), clusters.get(i));
					stDev2 = GetStDev(clusters.get(i), clusters.get(i + 1));
					stDevDiff1 = Math.abs(clusters.get(i - 1).getStDev() - clusters.get(i).getStDev());
					stDevDiff2 = Math.abs(clusters.get(i).getStDev() - clusters.get(i + 1).getStDev());
					countProp1 = (double)clusters.get(i - 1).getCount() / (double)clusters.get(i).getCount();
					if (countProp1 > 1.0)
					{
						countProp1 = 1.0 / countProp1;
					}
					countProp2 = (double)clusters.get(i).getCount() / (double)clusters.get(i + 1).getCount();
					if (countProp2 > 1.0)
					{
						countProp2 = 1.0 / countProp2;
					}

					if (diff1 < diff2 && diff1 <= maxDiff && stDev1 < maxStDev && stDevDiff1 < maxStDevDiff * countProp1)
					{
						joined.add(i);
						clusters.get(i - 1).Join(clusters.get(i));
						i = i + 2;
					}
					else if (diff1 >= diff2 && diff2 <= maxDiff && stDev2 < maxStDev && stDevDiff2 < maxStDevDiff * countProp2)
					{
						joined.add(i + 1);
						clusters.get(i).Join(clusters.get(i + 1));
						i = i + 3;
					}
					else
					{
						i = i + 1;
					}
				}
			}
			if (joined.isEmpty())
			{
				break;
			}
			for (int j = joined.size() - 1;j >= 0;j--)
			{
				clusters.remove(joined.get(j));
			}
		}
	}

	private void JoinInnerClusters(ArrayList<StatCluster> clusters)
	{
		int i;
		ArrayList<Integer> joined = new ArrayList<Integer>();
		while (true)
		{
			i = 1;
			joined.clear();
			while (i < clusters.size() - 1)
			{
				//if(!clusters[i-1].Close && !clusters[i-1].Close && !clusters[i-1].Close && clusters[i].Count < minCount && /*clusters[i-1].Count > minCount && clusters[i+1].Count > minCount && */clusters[i-1].StDev - clusters[i+1].StDev < maxDiff) {
				//if(!clusters[i-1].Close && !clusters[i-1].Close && !clusters[i-1].Close && clusters[i].Count < minCount && clusters[i-1].Count > minCount && clusters[i+1].Count > minCount && clusters[i-1].StDev - clusters[i+1].StDev < maxDiff) {
				if (!clusters.get(i - 1).getClose() && !clusters.get(i).getClose() && !clusters.get(i + 1).getClose() && clusters.get(i).getCount() < minCount && clusters.get(i - 1).getCount() + clusters.get(i + 1).getCount() > minCount && GetDiff(clusters.get(i - 1), clusters.get(i + 1)) < maxDiff && clusters.get(i).getStDev() < clusters.get(i).getStDev())
				{
					clusters.get(i - 1).Join(clusters.get(i));
					clusters.get(i - 1).Join(clusters.get(i + 1));
					joined.add(i);
					joined.add(i + 1);
					i = i + 3;
				}
				else
				{
					i = i + 1;
				}
			}
			if (joined.isEmpty())
			{
				break;
			}
			for (int j = joined.size() - 1;j >= 0;j--)
			{
				clusters.remove(joined.get(j));
			}
		}
	}

	private void CloseBigClusters(ArrayList<StatCluster> clusters)
	{
		for (StatCluster sc : clusters)
		{
			if (!sc.getClose() && sc.getCount() >= minCount)
			{
				sc.setClose(true);
			}
		}
	}

	//endregion

	//region Auxiliar Methods


	private double GetDiff(StatCluster sc1, StatCluster sc2)
	{
		return Math.abs(sc1.getMean() - sc2.getMean());
	}

	private double GetStDev(StatCluster sc1, StatCluster sc2)
	{
		double sum = sc1.getSum() + sc2.getSum();
		double sumSq = sc1.getSumSq() + sc2.getSumSq();
		double sqSum = sum * sum;
		int n = sc1.getCount() + sc2.getCount();
		double var = (sumSq - (sqSum / n)) / (n - 1);
		//double var = (n * sumSq - sum*sum) / (n * (n - 1));
		if (var < 1.0)
		{
			return 0.0;
		}
		return Math.sqrt(var);
	}

	//endregion

	//endregion

	//region Inner Class StatClust

	/**  Speciphical cluster class for Stationary Clustering 
	*/
	public static class StatCluster
	{

		//region Fields

		private int firstIndex;
		private int lastIndex;
		private int count;
		private double mean;
		private double biasMean;
		private double stDevMean;
		private double sum;
		private double sumSq;
		private boolean close;
		private double weight;
		private double minStDev;

		//endregion

		//region Constructors

		/**  Constructor 
		*/
		public StatCluster()
		{
		}

		/**  Constructor 
		 @param index index of cluster 
		 @param val value 
		*/
		public StatCluster(int index, double val)
		{
			firstIndex = index;
			lastIndex = index;
			count = 1;
			mean = val;
			sum = val;
			sumSq = val * val;
			weight = 0.0;
			minStDev = 0.5;
		}

		/**  Constructor 
		 @param sc cluster (copy constructor) 
		*/
		public StatCluster(StatCluster sc)
		{
			firstIndex = sc.getFirstIndex();
			lastIndex = sc.getLastIndex();
			count = sc.getCount();
			sum = sc.getSum();
			sumSq = sc.sumSq;
			weight = sc.getWeight();
			minStDev = 0.5;
		}

		/**  Constructor (join) 
		 @param sc1 first cluster 
		 @param sc2 second cluster 
		*/
		public StatCluster(StatCluster sc1, StatCluster sc2)
		{
			firstIndex = Math.min(sc1.getFirstIndex(), sc2.getFirstIndex());
			lastIndex = Math.max(sc1.getLastIndex(), sc2.getLastIndex());
			count = sc1.getCount() + sc2.getCount();
			sum = sc1.getSum() + sc2.getSum();
			sumSq = sc1.sumSq + sc2.getSumSq();
			weight = 0.0;
			minStDev = 0.5;
		}

		//endregion

		//region Properties

		/**  first index in the cluster 
		*/
		public final int getFirstIndex()
		{
			return firstIndex;
		}

		/**  last index in the cluster 
		*/
		public final int getLastIndex()
		{
			return lastIndex;
		}

		/**  number of elements 
		*/
		public final int getCount()
		{
			return count;
		}

		/**  mean of the cluster 
		*/
		public final double getMean()
		{
			return sum / (double)count;
		}

		/**  standard deviation of the cluster 
		*/
		public final double getStDev()
		{
			if (count == 1)
			{
				return minStDev;
			}
			else
			{
				double num = count * sumSq - sum * sum;
				if (num < 0.01)
				{
					return minStDev;
				}
				return Math.sqrt(num / (count * (count - 1)));
			}
		}

		/**  mean of the data bias 
		*/
		public final double getBiasMean()
		{
			return biasMean;
		}
		public final void setBiasMean(double value)
		{
			biasMean = value;
		}

		/**  standard deviation of the data bias 
		*/
		public final double getBiasStDev()
		{
			return stDevMean;
		}
		public final void setBiasStDev(double value)
		{
			stDevMean = value;
		}

		/**  precalculated sum of data (for standard deviation quick calculation) 
		*/
		public final double getSum()
		{
			return sum;
		}

		/**  precalculated squared sum of data (for standard deviation quick calculation)  
		*/
		public final double getSumSq()
		{
			return sumSq;
		}

		/**  if the cluster is closed 
		*/
		public final boolean getClose()
		{
			return close;
		}
		public final void setClose(boolean value)
		{
			close = value;
		}

		/**  cluster weight to be applied 
		*/
		public final double getWeight()
		{
			return weight;
		}
		public final void setWeight(double value)
		{
			weight = value;
		}

		//endregion

		//region Public Methods

		/**  Join this cluster with another 
		 @param sc the other cluster 
		*/
		public final void Join(StatCluster sc)
		{
			firstIndex = Math.min(firstIndex, sc.getFirstIndex());
			lastIndex = Math.max(lastIndex, sc.getLastIndex());
			count = count + sc.getCount();
			sum = sum + sc.getSum();
			sumSq = sumSq + sc.getSumSq();
			weight = weight + sc.weight;
		}

		//endregion

	}

	//endregion
}