package NumCalc;

//region Imports

import ClassicalStat.*;

//endregion


/**  Beta-Binomial Bayesian Model 
*/
public class BBModel  {

	//region Fields

	private BetaDist beta;
	private Statistics stat;
	private Combinatory comb;
	private RandomGen rand;

	//endregion

	//region Constructor

	/**  Constructor 
	*/
	public BBModel()
	{
		beta = new BetaDist();
		stat = new Statistics();
		comb = new Combinatory();
		rand = new RandomGen();
	}

	//endregion

	//region Properties

	/**  beta a parameter 
	*/
	public final double getA() 	{
		return beta.getA();
	}

	/**  beta b parameter 
	*/
	public final double getB() 	{
		return beta.getB();
	}


	//endregion

	//region Public Methods

	/**  Set informative prior with a and b defined 
	 @param a a parameter 
	 @param b b parameter 
	*/
	public final void SetInformativePrior(int a, int b)
	{
		beta.setA(a);
		beta.setB(b);
	}

	/**  Set uninformative prior Beta(1,1) equivalent to Uniform[0,1] 
	*/
	public final void SetUninformativePrior()
	{
		beta.setA(1);
		beta.setB(1);
	}

	public final void SetElicitedPrior(double m, double max)
	{
		int a;
		while (true)
		{
			a = rand.NextInt(0, 100);
			beta.setA(a);
			beta.SetB((int)beta.getA(), m);
			double p = beta.Probability(max, false);
			if (p > 0.95)
			{
				return;
			}
		}
	}

	/**  Bayesian conjugate update : Beta(a+x, b+n-x) 
	 @param x number of successes 
	 @param n number of trials 
	*/
	public final void Update(int x, int n)
	{
		beta.setA(beta.getA() + x);
		beta.setB(beta.getB() + n - x);
	}

	/**  Mean of proportion 
	 @return  the mean 
	*/
	public final double PropMean()
	{
		return beta.Mean();
	}

	/**  Variance of proportion 
	 @return  the variance 
	*/
	public final double PropVar()
	{
		return beta.Var();
	}


	/**  Predictive probablity of y successes in m trials while already x successes in n trials 
	 @param x new x successes 
	 @param n new n trials 
	 @param acum if it is acumulated probability or not 
	 @return  the probability 
	*/
	public final double Probability(int x, int n, boolean acum)
	{
		if (acum)
		{
			return ProbabilityAcum(x, n);
		}
		else
		{
			return Probability(x, n);
		}
	}

	private double Probability(int x, int n)
	{
		int a = (int)beta.getA();
		int b = (int)beta.getB();
		return (comb.Combinations(n, x) * stat.B(a + x, b + n - x)) / stat.B(a, b);
	}

	private double ProbabilityAcum(int x, int n)
	{
		double acum = 0;
		for (int i = 1; i <= x; i++)
		{
			acum += Probability(i, n);
		}
		return acum;
	}

	/**  Predictive probablity of y successes in m trials while already x successes in n trials 
	 @param y new y successes 
	 @param m new m trials 
	 @param x old x successes 
	 @param n old n trials 
	 @return 
	*/
	public final double Probability(int y, int m, int x, int n)
	{
		int a = (int) beta.getA();
		int b = (int) beta.getB();
		return (comb.Combinations(m, y) * stat.B(a + x + y, b + m + n - x - y)) / stat.B(a + x, b + n - x);
	}

	//endregion

	//region Private Methods


	//endregion
}