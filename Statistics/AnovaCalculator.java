package ClassicalStat;

import java.util.*;

public class AnovaCalculator
{

	//region Fields

	private ArrayList<Integer> groupCounts;

	private Statistics stat;
	private Combinatory comb;
	private ArrayList<Double> values;

	private ArrayList<Double> groupMeans;
	private ArrayList<Double> groupSs;

	private double ssInter;
	private double meanInter;
	private int dfInter;

	private double ssIntra;
	private double meanIntra;
	private int dfIntra;

	private double ssTotal;
	private double meanTotal;
	private int dfTotal;

	private double msInter;
	private double msIntra;
	private double F;

	private double pValue;
	private double significance;

	//endregion

	//region Constructor

	public AnovaCalculator()
	{
		stat = new Statistics();
		comb = new Combinatory();
	}

	//endregion

	//region Properties

	public final ArrayList<Integer> getGroupCounts()
	{
		return groupCounts;
	}
	public final void setGroupCounts(ArrayList<Integer> value)
	{
		groupCounts = value;
	}

	public final double getMeanTotal()
	{
		return meanTotal;
	}
	public final void setMeanTotal(double value)
	{
		meanTotal = value;
	}

	public final ArrayList<Double> getValues()
	{
		return values;
	}
	public final void setValues(ArrayList<Double> value)
	{
		values = value;
	}

	public final double getSSInter()
	{
		return ssInter;
	}

	public final double getSSIntra()
	{
		return ssIntra;
	}

	public final double getSSTotal()
	{
		return ssTotal;
	}

	public final int getDfInter()
	{
		return dfInter;
	}

	public final int getDfIntra()
	{
		return dfIntra;
	}

	public final double getMSInter()
	{
		return msInter;
	}

	public final double getMSIntra()
	{
		return msIntra;
	}

	public final double getFVal()
	{
		return F;
	}

	public final int getDfTotal()
	{
		return dfTotal;
	}

	public final double getPValue()
	{
		return pValue;
	}

	public final boolean getTestResult()
	{
		return pValue < significance;
	}

	//endregion

	//region Public Methods

	public final void Calculate()
	{
		 List<List<Double>> pGroups = GetPermutation(groupCounts, values);
		 CalcF(pGroups, values, meanTotal);
	}

	public final List<List<Double>> GetPermutation(ArrayList<Integer> groupCounts, ArrayList<Double> values)
	{
		ArrayList<Integer> valueIndexes = comb.PermutationNoRepeat(values.size());
		List<List<Double>> pGroups = new ArrayList<List<Double>>();
		List<Double> pGroup;
		int k = 0;
		for (int i = 0; i < groupCounts.size(); i++)
		{
			pGroup = new ArrayList<Double>();
			for (int j = 0; j < groupCounts.get(i); j++)
			{
				pGroup.add(values.get(valueIndexes.get(k++)));
			}
			pGroups.add(pGroup);
		}
		return pGroups;
	}

	public final double CalcF(List<List<Double>> groups, ArrayList<Double> values, double meanTotal)
	{
		this.values = values;
		dfTotal = values.size() - 1;
		dfInter = groups.size() - 1;
		dfIntra = dfTotal - dfInter;
		ArrayList<Double> groupMeans = new ArrayList<Double>();
		ArrayList<Double> groupSs = new ArrayList<Double>();
		for (List<Double> group : groups)
		{
			groupMeans.add(stat.Mean(group));
			groupSs.add(stat.SumSq(group));
		}

		ssInter = SSE(groupMeans, meanTotal, values.size() / groups.size());
		ssIntra = SST(groups, groupMeans);
		ssTotal = ssInter + ssIntra;

		if (ssInter == 0)
		{
			F = 0;
			return F;
		}
		if (ssIntra == 0)
		{
			F = Double.MAX_VALUE;
			return F;
		}

		msInter = ssInter / dfInter;
		msIntra = ssIntra / getDfIntra();

		F = msInter / msIntra;
		return F;
	}

	@Override
	public String toString()
	{
		StringBuilder sb = new StringBuilder();
		sb.append("Source\tSS\tdf\tMS\tF\tp");
		sb.append(ssInter + "\t" + getDfInter() + "\t" + getMSInter() + "\t" + F + "\t" + pValue);
		sb.append(ssIntra + "\t" + getDfIntra() + "\t" + getMSIntra());
		sb.append(ssTotal + "\t" + getDfTotal());
		return sb.toString();
	}

	//endregion

	//region Private Methods

	private double SSE(ArrayList<Double> groupMeans, double meanTotal, double groupCount)
	{
		double sse = 0;
		for (double gMean : groupMeans)
		{
			sse += Math.pow((gMean - meanTotal), 2);
		}
		return sse * groupCount;
	}

	private double SST(List<List<Double>> groups, ArrayList<Double> groupMeans)
	{
		double sst = 0;
		for (int i = 0; i < groups.size(); i++)
		{
			for (double val : groups.get(i))
			{
				sst = sst + Math.pow((val - groupMeans.get(i)), 2);
			}
		}
		return sst;
	}

	//endregion

}