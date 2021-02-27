package ClassicalStat;

//region Imports

import java.util.*;

//endregion


public class ANOVA {

	//region Fields

	private Statistics stat;
	private List<List<Double>> groups;
	private ArrayList<Double> data;

	private ArrayList<Double> groupMeans;
	private ArrayList<Double> groupSs;

	private double ssInter;
	private double meanInter;

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

	public ANOVA()
	{
		stat = new Statistics();
		data = new ArrayList<Double>();
		significance = 0.05;
	}

	//endregion

	//region Public Methods

	public final double Calculate(List<List<Double>> groups)
	{
		if (groups.size() <= 1)
		{
			pValue = 1;
			return pValue;
		}
		ArrayList<Double> values = new ArrayList<Double>();
		for (int i = 0; i < groups.size(); i++)
		{
			values.addAll(groups.get(i));
		}
		meanTotal = stat.Mean(values);

		AnovaCalculator ac = new AnovaCalculator();
		double F = ac.CalcF(groups, values, meanTotal);
		pValue = 1 - stat.FDist(F, ac.getDfInter(), ac.getDfIntra());
		return pValue;
	}

	public final double CalculateRand(List<List<Double>> groups, int iter, int maxThreads)
	{
		if (maxThreads <= 1)
		{
			return CalculateRandMono(groups, iter);
		}
		else
		{
			return -1;
		}
	}

	//endregion

	//region Private Methods

	private double CalculateRandMono(List<List<Double>> groups, int iter)
	{
		if (groups.size() <= 1)
		{
			pValue = 1;
			return pValue;
		}
		ArrayList<Integer> groupCounts = new ArrayList<Integer>();
		ArrayList<Double> values = new ArrayList<Double>();
		for (int i = 0; i < groups.size(); i++)
		{
			groupCounts.add(groups.get(i).size());
			values.addAll(groups.get(i));
		}
		meanTotal = stat.Mean(values);

		AnovaCalculator ac = new AnovaCalculator();
		double Fobs = ac.CalcF(groups, values, meanTotal);
		List<List<Double>> pGroups;
		double FPerm;
		ArrayList<Double> FPerms = new ArrayList<Double>();
		for (int i = 0; i < iter; i++)
		{
			pGroups = ac.GetPermutation(groupCounts, values);
			FPerm = ac.CalcF(pGroups, values, meanTotal);
			FPerms.add((double) (FPerm > Fobs ? 1 : 0));
		}
		pValue = stat.Mean(FPerms);
		return pValue;
	}

	//endregion

}