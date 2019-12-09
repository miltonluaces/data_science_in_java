package Clustering;

import Clustering.Properties.*;
import ClassicalStat.*;
import java.util.*;

//region Imports


//endregion


/**  Algorithmic repository class for stationarity time series clustering 
*/
public class HomogeneityClustering
{

	//region Fields

	private double origMaxDiff;
	private double origMaxStDev;
	private double maxDiff;
	private double maxStDev;
	private double maxStDevDiff;
	private int minCount;
	private HypoTest hypoTest;
	private double significance;
	private TestType testType = TestType.values()[0];
	private int minTestCount;

	//endregion

	//region Constructor

	/** 
	 Constructor 
	 
	 @param maxDiff maximum difference allowed 
	 @param maxStDev maximum standard deviation allowed 
	 @param minCount minimum number of elements allowed 
	 @param significance significance of homogeneity hypothesis test 
	 @param testType type of hypothesis test 
	*/
	public HomogeneityClustering(double maxDiff, double maxStDev, int minCount, double significance, TestType testType)
	{
		this.maxDiff = maxDiff;
		this.maxStDev = maxStDev;
		this.origMaxDiff = maxDiff;
		this.origMaxStDev = maxStDev;
		this.minCount = minCount;
		this.minTestCount = 5;
		this.hypoTest = new HypoTest(significance, minTestCount);
		this.significance = significance;
		this.testType = testType;
	}

	//endregion

	//region Properties

	public final double getSignificance()
	{
		return significance;
	}
	public final void setSignificance(double value)
	{
		significance = value;
	}

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

	/**  type of statistical test for joins 
	*/
	public final TestType getTest()
	{
		return testType;
	}
	public final void setTest(TestType value)
	{
		testType = value;
	}

	/**  cut-off for trimming 
	*/
	public final double getCutOff()
	{
		return hypoTest.CutOff;
	}
	public final void setCutOff(double value)
	{
		hypoTest.CutOff = value;
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
	public final ArrayList<StatCluster> Clustering(ArrayList<Double> serie, int firstIndex, int lastIndex)
	{
		if (serie == null)
		{
			throw new RuntimeException(Strings.Serie_is_null);
		}

		ArrayList<StatCluster> clusters = CreateClusters(serie, firstIndex, lastIndex);
		maxDiff = origMaxDiff;
		maxStDev = origMaxStDev;
		while (maxDiff < 100.0)
		{
			if (testType == TestType.Statistics)
			{
				JoinNeighbourClustersStatistics(clusters);
			}
			else
			{
				JoinNeighbourClusters(clusters);
			}
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

	/**  Test of homogeneity of samples 
	 @param serie time series 
	 @return  homogeneity of partitions 
	*/
	public final ArrayList<Double> TestHomogeneity(ArrayList<Double> serie)
	{
		testType = TestType.Wilcoxon;
		ArrayList<Double> homog = new ArrayList<Double>();
		StatCluster ini = new StatCluster();
		StatCluster end = new StatCluster();
		end.getValues().addAll(serie);
		for (int i = 0;i < minTestCount;i++)
		{
			ini.getValues().add(serie.get(i));
			end.getValues().remove(0);
			homog.add(0.0);
		}
		for (int i = minTestCount;i < serie.size() - minTestCount;i++)
		{
			homog.add(TestHomogeneity(ini, end));
			ini.getValues().add(serie.get(i));
			end.getValues().remove(0);
		}
		for (int i = 0;i < minTestCount;i++)
		{
			homog.add(0.0);
		}
		int changeModelIndex = -1;
		for (int i = homog.size() - minTestCount;i >= minTestCount;i--)
		{
			if (i > 0.5)
			{
				changeModelIndex = i;
			}
			break;
		}
		return homog;
	}

	/**  Test of homogeneity based in clusters 
	 @param clusters clusters to test 
	 @return  index of partition 
	*/
	public final int TestFirstHomogeneousCluster(ArrayList<StatCluster> clusters)
	{
		StatCluster ini;
		StatCluster end;
		for (int c = 0;c < clusters.size() - 1;c++)
		{
			ini = new StatCluster();
			end = new StatCluster();
			ini.getValues().add(clusters.get(c).getMean());
			for (int i = c + 1;i < clusters.size();i++)
			{
				end.getValues().add(clusters.get(i).getMean());
			}
			if (ini.getFirstIndex() > minTestCount && AreHomogeneous(ini, end))
			{
				return ini.getFirstIndex();
			}
		}
		return clusters.get(clusters.size() - 1).getFirstIndex();
	}

	//endregion

	//region Private Methods

	//region Calculate Methods

	private ArrayList<StatCluster> CreateClusters(ArrayList<Double> serie, int firstIndex, int lastIndex)
	{
		ArrayList<StatCluster> clusters = new ArrayList<StatCluster>();
		StatCluster sC;
		for (int j = firstIndex; j <= lastIndex; j++)
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
					if (AreHomogeneous(clusters.get(i - 1), clusters.get(i)))
					{
						joined.add(i);
						clusters.get(i - 1).Join(clusters.get(i));
						i = i + 2;
					}
					else if (AreHomogeneous(clusters.get(i), clusters.get(i + 1)))
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

	private void JoinNeighbourClustersStatistics(ArrayList<StatCluster> clusters)
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

	/**  Checks if two clusters are homogeneous 
	 @param sc1 cluster 1 
	 @param sc2 cluster 2 
	 @return  true if they are homogeneous, false if they are not 
	*/
	public final boolean AreHomogeneous(StatCluster sc1, StatCluster sc2)
	{
		double pValue = TestHomogeneity(sc1, sc2);
		return pValue >= significance;
	}

	/**  Checks if two clusters are homogeneous 
	 @param sc1 cluster 1 
	 @param sc2 cluster 2 
	 @param pValue return parameter pValue calculated 
	 @return  true if they are homogeneous, false if they are not 
	*/
	public final boolean AreHomogeneous(StatCluster sc1, StatCluster sc2, tangible.RefObject<Double> pValue)
	{
		pValue.argValue = TestHomogeneity(sc1, sc2);
		return pValue.argValue >= significance;
	}

	/**  Tests if two subsets divided by the index are homogeneous 
	 @param serie time series 
	 @param index index of first value of final subset 
	 @return  if the two subsets are homogeneous or not 
	*/
	public final boolean AreHomogeneous(ArrayList<Double> serie, int index)
	{
		double pValue = -1;
		tangible.RefObject<Double> tempRef_pValue = new tangible.RefObject<Double>(pValue);
		boolean tempVar = AreHomogeneous(serie, index, tempRef_pValue);
	pValue = tempRef_pValue.argValue;
	return tempVar;
	}

	/**  Tests if two subsets divided by the index are homogeneous 
	 @param serie time series 
	 @param index index of first value of final subset 
	 @param pValue return parameter of pValue calculated 
	 @return  if the two subsets are homogeneous or not 
	*/
	public final boolean AreHomogeneous(ArrayList<Double> serie, int index, tangible.RefObject<Double> pValue)
	{
		if (index <= 0 || index >= serie.size())
		{
			throw new RuntimeException(Strings.Index_must_be_positive_and_smaller_than_series_count);
		}

		StatCluster sc1 = new StatCluster();
		for (int i = 0; i < index; i++)
		{
			sc1.AddValueAndSum(serie.get(i));
		}
		StatCluster sc2 = new StatCluster();
		for (int i = index; i < serie.size(); i++)
		{
			sc2.AddValueAndSum(serie.get(i));
		}
		return AreHomogeneous(sc1, sc2, pValue);
	}


	/**  Tests if two subsets divided by the index are homogeneous 
	 @param serie time series 
	 @param iniIndex1 index of first value of initial subset 
	 @param endIndex1 index of last value of initial subset 
	 @param iniIndex2 index of first value of final subset 
	 @param endIndex2 index of last value of final subset 
	 @return  if the two subsets are homogeneous or not 
	*/
	public final boolean AreHomogeneous(ArrayList<Double> serie, int iniIndex1, int endIndex1, int iniIndex2, int endIndex2)
	{
		double pValue = -1;
		tangible.RefObject<Double> tempRef_pValue = new tangible.RefObject<Double>(pValue);
		boolean tempVar = AreHomogeneous(serie, iniIndex1, endIndex1, iniIndex2, endIndex2, tempRef_pValue);
	pValue = tempRef_pValue.argValue;
	return tempVar;
	}

	/**  Tests if two subsets divided by the index are homogeneous 
	 @param serie time series 
	 @param iniIndex1 index of first value of initial subset 
	 @param endIndex1 index of last value of initial subset 
	 @param iniIndex2 index of first value of final subset 
	 @param endIndex2 index of last value of final subset 
	 @param pValue return value of p value calcualted 
	 @return  if the two subsets are homogeneous or not 
	*/
	public final boolean AreHomogeneous(ArrayList<Double> serie, int iniIndex1, int endIndex1, int iniIndex2, int endIndex2, tangible.RefObject<Double> pValue)
	{
		if (iniIndex1 < 0 || iniIndex1 >= serie.size() || endIndex1 <= iniIndex1 || endIndex1 >= serie.size() || iniIndex2 <= iniIndex1 || iniIndex2 >= serie.size() || endIndex2 <= iniIndex2 || endIndex2 >= serie.size())
		{
			throw new RuntimeException(Strings.Bad_indexes);
		}

		StatCluster sc1 = new StatCluster();
		for (int i = iniIndex1;i < endIndex1;i++)
		{
			sc1.AddValueAndSum(serie.get(i));
		}
		StatCluster sc2 = new StatCluster();
		for (int i = iniIndex2;i < endIndex2;i++)
		{
			sc2.AddValueAndSum(serie.get(i));
		}
		return AreHomogeneous(sc1, sc2, pValue);
	}


	public final double TestHomogeneity(StatCluster sc1, StatCluster sc2)
	{
		if (sc1.getMean() == sc2.getMean() && sc1.getStDev() == sc2.getStDev())
		{
			return 1.0;
		}

		double pValue = -1;
		double pValueMean = -1;
		double pValueVar = -1;
		switch (testType)
		{
			case ParametricMean:
				tangible.RefObject<Double> tempRef_pValueMean = new tangible.RefObject<Double>(pValueMean);
				tangible.RefObject<Double> tempRef_pValueVar = new tangible.RefObject<Double>(pValueVar);
				hypoTest.ParametricHomogeneity(sc1.getValues(), sc2.getValues(), tempRef_pValueMean, tempRef_pValueVar);
			pValueVar = tempRef_pValueVar.argValue;
			pValueMean = tempRef_pValueMean.argValue;
				pValue = pValueMean;
				break;
			case ParametricVar:
				tangible.RefObject<Double> tempRef_pValueMean2 = new tangible.RefObject<Double>(pValueMean);
				tangible.RefObject<Double> tempRef_pValueVar2 = new tangible.RefObject<Double>(pValueVar);
				hypoTest.ParametricHomogeneity(sc1.getValues(), sc2.getValues(), tempRef_pValueMean2, tempRef_pValueVar2);
			pValueVar = tempRef_pValueVar2.argValue;
			pValueMean = tempRef_pValueMean2.argValue;
				pValue = pValueVar;
				break;
			case ParametricMeanVar:
				tangible.RefObject<Double> tempRef_pValueMean3 = new tangible.RefObject<Double>(pValueMean);
				tangible.RefObject<Double> tempRef_pValueVar3 = new tangible.RefObject<Double>(pValueVar);
				hypoTest.ParametricHomogeneity(sc1.getValues(), sc2.getValues(), tempRef_pValueMean3, tempRef_pValueVar3);
			pValueVar = tempRef_pValueVar3.argValue;
			pValueMean = tempRef_pValueMean3.argValue;
				pValue = Math.max(pValueMean, pValueVar);
				break;
			case MannWithney:
				pValue = hypoTest.MannWhitney(sc1.getValues(), sc2.getValues(), true);
				break;
			case Wilcoxon:
				pValue = hypoTest.WaldWolfowitz(sc1.getValues(), sc2.getValues(), true);
				break;
		}
		return pValue;
	}

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

	//region Inner Class StatCluster

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
		private ArrayList<Double> values;

		//endregion

		//region Constructors

		/**  Constructor 
		*/
		public StatCluster()
		{
			values = new ArrayList<Double>();
			sum = 0;
			sumSq = 0;
			weight = 1;
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
			values = new ArrayList<Double>();
			values.add(val);
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
			values = new ArrayList<Double>();
			values.addAll(sc.getValues());
			sum = 0;
			sumSq = 0;
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
			values = new ArrayList<Double>();
			values.addAll(sc1.getValues());
			values.addAll(sc2.getValues());
		}

		//endregion

		//region Properties

		/**  first index in the cluster 
		*/
		public final int getFirstIndex()
		{
			return firstIndex;
		}
		public final void setFirstIndex(int value)
		{
			firstIndex = value;
		}

		/**  last index in the cluster 
		*/
		public final int getLastIndex()
		{
			return lastIndex;
		}
		public final void setLastIndex(int value)
		{
			lastIndex = value;
		}

		/**  number of elements 
		*/
		public final int getCount()
		{
			return count;
		}
		public final void setCount(int value)
		{
			count = value;
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

		/**  list of cluster values 
		*/
		public final ArrayList<Double> getValues()
		{
			return values;
		}
		public final void setValues(ArrayList<Double> value)
		{
			values = value;
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
			count = lastIndex - firstIndex + 1;
			sum = 0;
			sumSq = 0;
			for (double val : sc.values)
			{
				sum += val;
				sumSq += val * val;

			}
			for (double val : this.values)
			{
				sum += val;
				sumSq += val * val;
			}
			weight = weight + sc.weight;
			values.addAll(sc.getValues());
		}

		/**  Add one value to the cluster 
		 @param val the value 
		*/
		public final void AddValue(double val)
		{
			this.values.add(val);
			this.count = this.count + 1;
		}

		/**  Add one value to the cluster and increment sums 
		 @param val the value 
		*/
		public final void AddValueAndSum(double val)
		{
			this.values.add(val);
			this.sum = this.sum + val;
			this.sumSq = this.sumSq + val * val;
			this.count = this.count + 1;
		}

		public final void Update(ArrayList<Double> newValues)
		{
			this.values = newValues;
			this.count = values.size();
			this.sum = 0;
			this.sumSq = 0;
			for (double val : newValues)
			{
				sum += val;
				this.sumSq = this.sumSq + val * val;
			}
		}

		//endregion

	}

	//endregion

	//region Enum TestType

	/**  Type of hypothesis testing for homogeneity 
	*/
	public enum TestType
	{
		/**  Join by statistics 
		*/
		Statistics,
		/**  Parametric test for mean 
		*/
		ParametricMean,
		/**  Parametric test for var 
		*/
		ParametricVar,
		/**  Parametric test for both mean and var 
		*/
		ParametricMeanVar,
		/**  Mann-Withney non parametric test 
		*/
		MannWithney,
		/**  Wilcoxon non parametric test 
		*/
		Wilcoxon;

		public int getValue()
		{
			return this.ordinal();
		}

		public static TestType forValue(int value)
		{
			return values()[value];
		}
	}

	//endregion
}