package NumCalc;

import ClassicalStat.*;
import java.util.*;

//region Imports


//endregion


/**  Dirichlet-Multinomial Bayesian Model (Polya distribution model) 
*/
public class DMModel
{

	//region Fields

	private DirichletDist dirichlet;
	private HashMap<String, Option> dirichletByName;
	private HashMap<Integer, Option> dirichletByIndex;
	private Statistics stat;
	private Combinatory comb;
	private RandomGen rGen;
	private int maxIt;

	//endregion

	//region Constructor

	/**  Constructor 
	*/
	public DMModel(int maxIt)
	{
		this.dirichlet = new DirichletDist();
		this.dirichletByName = new HashMap<String, Option>();
		this.dirichletByIndex = new HashMap<Integer, Option>();
		this.stat = new Statistics();
		this.comb = new Combinatory();
		this.rGen = new RandomGen();
		this.maxIt = maxIt;
	}

	//endregion

	//region Properties

	/**  Discount Factor 
	*/
	public final double getDiscFactor()
	{
		return dirichlet.getDiscFactor();
	}
	public final void setDiscFactor(double value)
	{
		dirichlet.setDiscFactor(value);
	}

	/**  Standard update quantity 
	*/
	public final double getUpdatePeriod()
	{
		return dirichlet.getUpdatePeriod();
	}
	public final void setUpdatePeriod(double value)
	{
		dirichlet.setUpdatePeriod(value);
	}

	/**  Real total quantity entered 
	*/
	public final double getRealTotalQty()
	{
		return dirichlet.getRealTotA();
	}
	public final void setRealTotalQty(double value)
	{
		dirichlet.setRealTotA(value);
	}

	//endregion

	//region Setters and Getters

	// region Options

	/**  Add a new option to the model 
	 @param name options name 
	*/
	public final void AddOption(String name)
	{
		int index = dirichletByName.size();
		Option opt = new Option(index, name, -1);
		dirichletByName.put(name, opt);
		dirichletByIndex.put(index, opt);
	}

	public final int GetOptionIndex(String name)
	{
		if (!dirichletByName.containsKey(name))
		{
			throw new RuntimeException("Error. Option " + name + " does not exist");
		}
		return dirichletByName.get(name).index;
	}

	public final void Mask(String name, boolean mask)
	{
		if (!dirichletByName.containsKey(name))
		{
			throw new RuntimeException("Error. Option " + name + " does not exist");
		}
		dirichletByName.get(name).mask = mask;
	}

	/**  Get a list of option names 
	 @return  the list of names 
	*/
	public final ArrayList<String> GetNames()
	{
		ArrayList<String> names = new ArrayList<String>();
		for (String name : dirichletByName.keySet())
		{
			names.add(name);
		}
		return names;
	}

	/**  To String, showing parameter names and its proportions 
	 @return  the string to show 
	*/
	@Override
	public String toString()
	{
		String str = "";
		for (Option opt : dirichletByName.values())
		{
			str += opt.name + " : " + dirichlet.Mean(opt.index) + "\n";
		}
		return str;
	}

	/**  Get a proportion of a particular index 
	 @param index the index 
	 @return  the proportion 
	*/
	public final double GetProportion(int index)
	{
		return dirichlet.Mean(index);
	}

	/**  Get a proportion of a particular name 
	 @param name the name 
	 @return  the proportion 
	*/
	public final double GetProportion(String name)
	{
		int i = dirichletByName.get(name).index;
		return dirichlet.Mean(i);
	}

	/**  Get a list of proportions 
	 @return  the list 
	*/
	public final ArrayList<Double> GetProportions()
	{
		ArrayList<Double> proportions = new ArrayList<Double>();
		for (Option opt : dirichletByName.values())
		{
			proportions.add(dirichlet.Mean(opt.index));
		}
		return proportions;
	}

	public final int GetNParams()
	{
		return dirichletByName.size();
	}

	//endregion

	//region Priors

	/**  Set an uninformative prior for the dirichlet 
	*/
	public final void SetUnInformativePrior()
	{
		for (int i = 0; i < dirichletByName.values().size();i++)
		{
			dirichlet.AddA(-1);
		}
		for (Option opt : dirichletByName.values())
		{
			opt.a = 1.0;
			if (opt.mask)
			{
				dirichlet.SetA(opt.index, 0);
			}
			else
			{
				dirichlet.SetA(opt.index, 1.0);
			}
		}
	}

	/**  Set an informative prior for the dirichlet 
	 @param successes successes for each option 
	 @param n number of trials done (or estimated) for proportion calculations 
	*/
	public final void SetInformativePrior(ArrayList<Integer> successes)
	{
		if (successes.size() != dirichlet.GetA().size())
		{
			throw new RuntimeException("Error. Number of proportios must fit dirichlet distribution parameters");
		}
		for (int i = 0; i < dirichletByName.values().size(); i++)
		{
			dirichlet.AddA(-1);
		}
		for (Option opt : dirichletByName.values())
		{
			opt.a = successes.get(opt.index);
			if (opt.mask)
			{
				dirichlet.SetA(opt.index, 0);
			}
			else
			{
				dirichlet.SetA(opt.index, successes.get(opt.index));
			}
		}
	}

	//endregion

	//endregion

	//region Public Methods


	/**  Predictive distribution [ G(Sum a_i) / G(Sum a_i + Sum x_i) ] * Prod [ G(a_i + x_i) / G(a_i)] 
	 @param X X successes 
	 @param cdf if it is cumulated distribution funcition (if not its pdf) 
	 @return  the probability 
	*/
	public final double Probability(List<Integer> X, boolean cdf)
	{
		if (cdf)
		{
			return ProbabilityCdf(X);
		}
		else
		{
			return Probability(X);
		}
	}

	private double Probability(List<Integer> X)
	{
		//double sA = dirichlet.GetTotA();
		double sA = this.getRealTotalQty();
		double sX = 0;
		for (int i = 0;i < X.size();i++)
		{
			sX += X.get(i);
		}

		double term1 = comb.Factorial((int)sX) * stat.G(sA) / stat.G(sX + sA);

		double term2 = 1;
		for (int i = 0; i < X.size(); i++)
		{
			double ai = dirichlet.GetA(i);
			term2 = term2 * stat.G(ai + X.get(i)) / (stat.G(ai) * comb.Factorial(X.get(i)));
		}
		return term1 * term2;
	}


   private double ProbabilityCdf(List<Integer> X)
   {
		int totX = 0;
		for (int x : X)
		{
			totX += x;
		}
		ArrayList<Integer> Xr;
		ArrayList<Double> Pr = new ArrayList<Double>();
		for (int it = 0; it < maxIt; it++)
		{
			Xr = new ArrayList<Integer>();
			for (int i = 0; i < X.size(); i++)
			{
				if (dirichletByIndex.get(i).mask)
				{
					Xr.add(0);
				}
				else
				{
					Xr.add(rGen.NextInt(0, X.get(i)));
				}
			}
			Pr.add(Probability(Xr));
		}
		double p = stat.Mean(Pr);
		return p;
   }

	/**  Quantile of each option for a particular probability 
	 @param p probability 
	 @param totX total of all options 
	 @return  a vector with option quantiles 
	*/
	public final ArrayList<Integer> Quantile(DotNetHelpers.RefObject<Double> p, DotNetHelpers.RefObject<Integer> totX)
	{

		ArrayList<Trial> trials = new ArrayList<Trial>();
		double pr;
		ArrayList<Integer> Xr;
		Trial trial;
		int tot = totX.argValue;
		int x;
		double totPr = 0;
		for (int it = 0; it < maxIt; it++)
		{

			Xr = GetTrial(totX.argValue);
			pr = Probability(Xr);
			totPr += pr;
			trial = new Trial(Xr, pr);
			trials.add(trial);

		}
		Collections.sort(trials);
		ArrayList<Integer> max = SwappingByMax(trials, p, totX);
		return max;
	}

	/**  Bayesian update of the dirichlet distribution for time series 
	 @param XSeries time series for update (all with the same count, fill with zeros otherwise) 
	 @param period update period 
	*/
	public final void Update(HashMap<String, ArrayList<Integer>> XSeries, int period)
	{
		ArrayList<Integer> serie;
		int firstIndex = 0;
		int lastIndex = period;
		boolean end = false;
		while (!end)
		{
			if (lastIndex >= XSeries.size() - 1)
			{
				lastIndex = XSeries.size() - 1;
				end = true;
			}
			HashMap<String, Integer> XDict = new HashMap<String, Integer>();
			double total;
			int totalInt;
			for (String optName : XSeries.keySet())
			{
				serie = XSeries.get(optName);
				total = 0;
				for (int i = firstIndex; i <= lastIndex; i++)
				{
					total += serie.get(i);
				}
				totalInt = (int)(Math.round(total));
				XDict.put(optName, totalInt);
			}
			Update(XDict, lastIndex - firstIndex);
			firstIndex = lastIndex + 1;
			lastIndex = firstIndex + 1 + period;
		}
	}

	/**  Bayesian update of the dirichlet distribution 
	 @param X values for update 
	 @param period update period 
	*/
	public final void Update(HashMap<String, Integer> XDict, int period)
	{
		//calculate update factor
		double n = period / dirichlet.getUpdatePeriod();
		double factor = Math.pow((1.00 - dirichlet.getDiscFactor()), n);

		//calculate new real total of A and normaliza A
		double totalQty = 0;
		for (int qty : XDict.values())
		{
			totalQty += qty;
		}
		dirichlet.setRealTotA(dirichlet.getRealTotA() * factor);
		dirichlet.NormalizeA();

		//normalize dictionary of new data           
		double strFactor = (totalQty * dirichlet.getMaxTotA()) / dirichlet.getRealTotA();
		HashMap<String, Integer> NormXDict = new HashMap<String, Integer>();
		for (String optName : XDict.keySet())
		{
			NormXDict.put(optName, (int)(Math.round(XDict.get(optName) * strFactor)));
		}

		//update distribution
		Option opt;
		for (String optName : NormXDict.keySet())
		{
			opt = dirichletByName.get(optName);
			if (opt.mask)
			{
				continue;
			}
			dirichlet.UpdateA(opt.index, XDict.get(optName));
		}
	}

	/**  Bayesian update of the dirichlet distribution 
	 @param X values for update 
	 @param period update period 
	*/
	public final void Update(List<Integer> X, int period)
	{
		HashMap<String, Integer> XDict = new HashMap<String, Integer>();
		Option opt;
		for (int i = 0; i < X.size(); i++)
		{
			opt = dirichletByIndex.get(i);
			if (!opt.mask)
			{
				XDict.put(opt.name, X.get(i));
			}
		}
		Update(XDict, period);
	}

	/**  Mean of the theta_i parameter 
	 @param i index 
	 @return  the mean 
	*/
	public final double Mean(int i)
	{
		return dirichlet.Mean(i);
	}

	/**  Covariance of the theta_i parameter 
	 @param i index 1 
	 @param j index 2 
	 @return  the covariance 
	*/
	public final double Cov(int i, int j)
	{
		return dirichlet.Cov(i, j);
	}

	/**  Marginal probability for one option and a particular value of theta 
	 @param i option index 
	 @param theta value of the parameter 
	 @param acum if it is accumulated probability or not 
	 @return  the probability 
	*/
	public final double ProbabilityMarg(int i, double theta, boolean acum)
	{
		return dirichlet.ProbabilityMargBeta(i, theta, acum);
	}

	/**  Marginal quantile for one option and a particular probability 
	 @param i option index 
	 @param p probability 
	 @return  the quantile 
	*/
	public final double QuantileMarg(int i, double p)
	{
		return dirichlet.QuantileMargBeta(i, p);
	}

	//endregion

	//region Private Methods

	//region Quantile methods

	private ArrayList<Integer> GetTrial(int total)
	{
		int n = GetNParams();
		ArrayList<Integer> X = new ArrayList<Integer>();
		int x;
		for (int i = 0; i < n; i++)
		{
			if (dirichletByIndex.get(i).mask)
			{
				x = 0;
			}
			else
			{
				x = rGen.NextInt(0, total);
			}
			X.add(x);
			total -= x;
		}
		if (total > 0)
		{
			int index = rGen.NextInt(0, n - 1);
			X.set(index, X.get(index) + total);
		}
		ArrayList<Integer> perm = comb.Permutation(0,n - 1);
		ArrayList<Integer> XPerm = new ArrayList<Integer>();
		for (int i : perm)
		{
			XPerm.add(X.get(i));
		}
		for (int i = 0; i < n; i++)
		{
			if (dirichletByIndex.get(i).mask)
			{
				XPerm.set(i, 0);
			}
		}
		return XPerm;
	}

	private int GetIndex(ArrayList<Trial> trials, DotNetHelpers.RefObject<Double> p, DotNetHelpers.RefObject<Integer> qty)
	{
		double totProb = 0;
		for (int i = 0; i < trials.size();i++)
		{
			totProb += trials.get(i).p;
		}

		double prob = 0;

		double quantity = 0;
		for (int i = 0; i < trials.size(); i++)
		{
			prob += trials.get(i).p;
			quantity += trials.get(i).getTotal();
			double pr = prob / totProb;
			if (p.argValue > 0)
			{
				if (prob / totProb >= p.argValue)
				{
					return i;
				}
			}
			else
			{
				if (quantity >= qty.argValue)
				{
					return i;
				}
			}
		}
		return trials.size() - 1;
	}

	private ArrayList<Integer> GetMax(ArrayList<Trial> trials, int index)
	{
		int n = GetNParams();
		ArrayList<Integer> max = new ArrayList<Integer>();
		for (int i = 0; i < n; i++)
		{
			max.add(0);
		}

		Trial trial;
		for (int i = 0; i <= index; i++)
		{
			trial = trials.get(i);
			for (int j = 0; j < n; j++)
			{
				if (trial.x.get(j).compareTo(max.get(j)) > 0)
				{
					max.set(j, trial.x.get(j));
				}
			}
		}
		return max;
	}

	private void Swap(ArrayList<Trial> trials, int i, int j)
	{
		Trial aux = trials.get(i);
		trials.set(i, trials.get(j));
		trials.set(j, aux);
	}


	private int GetTrial(ArrayList<Trial> trials, int ini, int end, ArrayList<Integer> max, boolean lower)
	{
		for (int i = ini; i <= end; i++)
		{
			if (lower)
			{
				if (trials.get(i).IsLowerThan(max))
				{
					return i;
				}
			}
			else
			{
				if (trials.get(i).IsHigherThan(max))
				{
					return i;
				}
			}


		}
		return -1;
	}

	private boolean SwapByMax(ArrayList<Trial> trials, ArrayList<Integer> max, int maxIndex)
	{
		int lowerIndex = GetTrial(trials, maxIndex + 1, trials.size() - 1, max, true);
		if (lowerIndex == -1)
		{
			return false;
		}
		for (int i = maxIndex; i > 0; i--)
		{
			if (trials.get(i).IsHigherThan(max))
			{
				Swap(trials, i, lowerIndex);
				return true;
			}
		}
		return false;
	}

	private ArrayList<Integer> SwappingByMax(ArrayList<Trial> trials, DotNetHelpers.RefObject<Double> p, DotNetHelpers.RefObject<Integer> qty)
	{
		int maxIndex = GetIndex(trials, p, qty);
		ArrayList<Integer> max = GetMax(trials, maxIndex);
		while (true)
		{
			boolean changed = SwapByMax(trials, max, maxIndex);
			if (!changed)
			{
				return max;
			}
			maxIndex = GetIndex(trials, p, qty);
			max = GetMax(trials, maxIndex);
		}
	}

	//endregion

	@Deprecated
	private void NormalizeProportions(ArrayList<Double> proportions)
	{
		double total = 0;
		for (double p : proportions)
		{
			total += p;
		}
		if (total != 1)
		{
			double newTotal = 0;
			for (int i = 0; i < proportions.size();i++)
			{
				proportions.set(i, proportions.get(i) / total);
				newTotal += proportions.get(i);
				if (newTotal > 1.0)
				{
					proportions.set(i, proportions.get(i) - newTotal - 1.0);
				}
			}
		}
	}

	//endregion

	//region Class Option

	/**  class Option for each dimension of the dirichlet 
	*/
	public static class Option
	{

		/** @param index index of the dimension 
		*/
		public int index;
		/** @param name name of the dimension 
		*/
		public String name;
		/** @param a a value for the dimension 
		*/
		public double a;
		/**  if this option should be not considered  
		*/
		public boolean mask;
		/**  Constructor 
		 @param index index of the dimension 
		 @param name name of the dimension 
		 @param a a value for the dimension 
		*/
		public Option(int index, String name, double a)
		{
			this.index = index;
			this.name = name;
			this.a = a;
			this.mask = false;
		}
	}

	//endregion

	//region Class Trial 

	public static class Trial implements java.lang.Comparable
	{
		public ArrayList<Integer> x;
		public double p;

		public Trial(ArrayList<Integer> x, double p)
		{
			this.x = x;
			this.p = p;
		}

		public final double getTotal()
		{
			double total = 0;
			for (int val : x)
			{
				total += val;
			}
			return total;
		}

		public final boolean IsLowerThan(ArrayList<Integer> max)
		{
			for (int i = 0; i < x.size(); i++)
			{
				if (x.get(i).compareTo(max.get(i)) > 0)
				{
					return false;
				}
			}
			return true;
		}

		public final boolean IsHigherThan(ArrayList<Integer> max)
		{
			for (int i = 0; i < x.size(); i++)
			{
				if (x.get(i).compareTo(max.get(i)) <= 0)
				{
					return false;
				}
			}
			return true;
		}

		public final int compareTo(Object obj)
		{
			if (this.p > ((Trial)obj).p)
			{
				return -1;
			}
			else if (this.p < ((Trial)obj).p)
			{
				return 1;
			}
			else
			{
				return 0;
			}
		}
	}

	//endregion

}