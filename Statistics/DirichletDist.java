package ClassicalStat;

import java.util.*;

//region Imports


//endregion


/**  Dirichlet distribution 
*/
public class DirichletDist
{

	//region Fields

	//a_i parameters 
	private ArrayList<Double> A;
	private double totA;
	private double realTotA;
	private Statistics stat;
	private ArrayList<BetaDist> bd;
	private double discFactor;
	private double updatePeriod;
	private double maxTotA;

	//endregion

	//region Constructor

	/**  Constructor 
	*/
	public DirichletDist()
	{
		stat = new Statistics();
		bd = new ArrayList<BetaDist>();
		A = new ArrayList<Double>();
		totA = 0;
		discFactor = 0.0;
		updatePeriod = 1.0;
		maxTotA = 100;
	}

	//endregion

	//region Properties

	/**  real total A 
	*/
	public final double getRealTotA()
	{
		return realTotA;
	}
	public final void setRealTotA(double value)
	{
		realTotA = value;
	}

	/**  Discount factor (between 0 and 1) for update 
	*/
	public final double getDiscFactor()
	{
		return discFactor;
	}
	public final void setDiscFactor(double value)
	{
		discFactor = value;
	}

	/**  standard period (in days) for update 
	*/
	public final double getUpdatePeriod()
	{
		return updatePeriod;
	}
	public final void setUpdatePeriod(double value)
	{
		updatePeriod = value;
	}

	/**  maximum quantity of totA before normalization 
	*/
	public final double getMaxTotA()
	{
		return maxTotA;
	}
	public final void setMaxTotA(double value)
	{
		maxTotA = value;
	}

	//endregion

	//region Setters and Getters

	/**  Clear previous information 
	*/
	public final void Clear()
	{
		A.clear();
		totA = 0;
	}

	/**  add a new dimension for A parameter 
	*/
	public final int AddA(double a)
	{
		A.add(a);
		totA += a;
		if (totA > maxTotA)
		{
			NormalizeA();
		}

		BetaDist b = new BetaDist();
		b.setA(a);
		bd.add(b);
		SetMarginal(bd.size() - 1);

		return A.size() - 1;
	}

	/**  set an existing a parameter 
	*/
	public final void SetA(int i, double a)
	{
		totA -= A.get(i);
		A.set(i, a);
		totA += A.get(i);
		if (totA > maxTotA)
		{
			NormalizeA();
		}

		SetMarginal(i);
	}

	/**  get value of an existing a parameter 
	*/
	public final double GetA(int i)
	{
		return A.get(i);
	}

	/**  get total of A vector 
	*/
	public final double GetTotA()
	{
		return totA;
	}

	/**  bayesian update of each element of the A vector 
	*/
	public final void UpdateA(int i, double a)
	{
		A.set(i, A.get(i) + a);
		totA += a;

		SetMarginal(i);
	}

	/**  the A vector 
	 @return  A vector 
	*/
	public final ArrayList<Double> GetA()
	{
		return A;
	}

	//endregion

	//region Public Methods

	//region Statistics

	/**  Expectation E[theta_i] 
	 @param i index 
	 @return  mean 
	*/
	public final double Mean(int i)
	{
		return A.get(i) / totA;
	}

	/**  Covariance(i,j) = a_i a_i / (/sum ai)^2 (sum ai+1) 
	 @param i i index 
	 @param j j index 
	 @return  the covariance value 
	*/
	public final double Cov(int i, int j)
	{
		double sum = 0;
		for (double a : A)
		{
			sum += a;
		}
		double num = A.get(i) * A.get(j);
		double den = Math.pow(sum, 2) * (sum + A.size());
		return num / den;
	}

	/**  a and b parameters of Beta_i(a,b) 
	 @param i index 
	 <return> Beta distribution 
	*/
	public final BetaDist GetBeta(int i)
	{
		if (i > bd.size() - 1)
		{
			throw new RuntimeException("Error. There are only " + i + " marginal beta distributions.");
		}

		return bd.get(i);
	}

	/**  Probability for a marginal beta 
	 @param i index 
	 @param theta proportion 
	 @param accum if it is accumulated or not 
	 @return 
	*/
	public final double ProbabilityMargBeta(int i, double theta, boolean accum)
	{
		if (i > bd.size() - 1)
		{
			throw new RuntimeException("Error. There are only " + i + " marginal beta distributions.");
		}

		return bd.get(i).Probability(theta, accum);
	}

	/**  Marginal beta distribution for the dirichlet 
	 @param i index 
	 @param p proportion 
	 @return 
	*/
	public final double QuantileMargBeta(int i, double p)
	{
		if (i > bd.size() - 1)
		{
			throw new RuntimeException("Error. There are only " + i + " marginal beta distributions.");
		}

		return bd.get(i).Quantile(p);
	}

	//endregion

	//region Main Functions


	/**  Probability function : 1/B(A) . Prod(theta_i^{a_i-1}) 
	 @param Theta vectorial parameter Theta 
	 @return  the probability 
	*/
	public final double Probability(ArrayList<Double> Theta)
	{
		double prodThetaA = 1;
		for (int i = 0;i < Theta.size();i++)
		{
			prodThetaA *= Math.pow(Theta.get(i),A.get(i) - 1);
		}
		double num = prodThetaA;
		double den = stat.B(A);
		return num / den;
	}

	//endregion

	//region Private Methods

	private void SetMarginal(int i)
	{
		double b = 0;
		for (int j = 0; j < A.size();j++)
		{
			if (j != i)
			{
				b += A.get(j);
			}
		}

		if (i > bd.size() - 1)
		{
			throw new RuntimeException("Error. There are only " + i + " marginal beta distributions.");
		}
		bd.get(i).setB(b);
	}

	public final void NormalizeA()
	{
		double newTotA = 0;
		double newA;
		for (int i = 0; i < A.size(); i++)
		{
			newA = this.Mean(i) * maxTotA;
			A.set(i, newA);
			newTotA += newA;
		}
		totA = newTotA;
	}

	public final void ApplyDiscFactor(double period)
	{
		double n = period / updatePeriod;
		double factor = Math.pow((1.00 - discFactor), n);
		for (int i = 0; i < A.size(); i++)
		{
			A.set(i, A.get(i) * factor);
		}
		totA *= factor;
	}

	//endregion

	//endregion

}