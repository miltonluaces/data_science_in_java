package ClassicalStat;

/**  Class for calculations of probability, accumulated probabiliy and statistics of discrete distributions 
*/
public class DiscreteDist
{

	//region Fields

	private Combinatory comb;

	//endregion

	//region Constructors

	/**  Constructor 
	*/
	public DiscreteDist()
	{
		comb = new Combinatory();
	}

	//endregion

	//region Public Methods

	//region Probability

	/**  Binomial Distribution 
	 Bernouilli Distribution when n = 1
	 @param x number of successes 
	 @param n number of independent yes/no experiments 
	 @param p probability of a success 
	 @param acum accumulated distribution or not 
	 @return  probability of x successes 
	*/
	public final double Binomial(int x, int n, double p, boolean acum)
	{
		if (acum)
		{
			double sum = 0.0;
			for (int i = 0;i <= x;i++)
			{
				sum += Binomial(i, n, p);
			}
			return sum;
		}
		else
		{
			return Binomial(x, n, p);
		}
	}

	/** 
	 
	 
	 @param x number of successes 
	 @param N number of independent yes/no experiments of the population 
	 @param n number of independent yes/no experiments of the sample 
	 @param d number of favorable elements in the sample 
	 @param acum accumulated distribution or not 
	 @return  probability of x successes 
	*/
	public final double Hypergeometric(int x, int N, int n, int d, boolean acum)
	{
		if (acum)
		{
			double sum = 0.0;
			for (int i = 0;i <= x;i++)
			{
				sum += Hypergeometric(i, N, n, d);
			}
			return sum;
		}
		else
		{
			return Hypergeometric(x, N, n, d);
		}
	}

	/**  Poisson Distribution 
	 @param x number of occurrences of an event 
	 @param lambda the expected number of occurrences that occur during the given interval 
	 @param acum accumulated distribution or not 
	 @return  probability of x occurrences 
	*/
	public final double Poisson(int x, double lambda, boolean acum)
	{
		if (acum)
		{
			double sum = 0.0;
			for (int i = 0;i <= x;i++)
			{
				sum += Poisson(i, lambda);
			}
			return sum;
		}
		else
		{
			return Poisson(x, lambda);
		}
	}

	/**  Hurdle Poisson Distribution 
	 @param x number of occurrences of an event 
	 @param alpha the probability of no occurrences during the given interval 
	 @param lambda the expected number of occurrences that occur during the given interval 
	 @param acum accumulated distribution or not 
	 @return  probability of x occurrences 
	*/
	public final double HurdlePoisson(int x, double alpha, double lambda, boolean acum)
	{
		if (x > 160)
		{
			x = 160;
		} //max value
		if (acum)
		{
			double sum = 0.0;
			for (int i = 0; i <= x; i++)
			{
				sum += HurdlePoisson(i, alpha, lambda);
			}
			return sum;
		}
		else
		{
			return HurdlePoisson(x, alpha, lambda);
		}
	}



	/** 
	 
	 
	 @param x number of experiments (fails) needed until n is reached 
	 Geometric / Pascal / Poyla distribution when n = 1
	 @param n specified number of succeses 
	 @param p probability of a success 
	 @param acum accumulated distribution or not 
	 @return  probability of x fails needed 
	*/
	public final double NegativeBinomial(int x, int n, double p, boolean acum)
	{
		if (acum)
		{
			double sum = 0.0;
			for (int i = 0;i <= x;i++)
			{
				sum += NegativeBinomial(i, n, p);
			}
			return sum;
		}
		else
		{
			return NegativeBinomial(x, n, p);
		}
	}

	//endregion

	//region Quantil

	/**  Inverse Binomial Distribution 
	 @param prob accumulated probabiliy of successes 
	 @param n number of independent yes/no experiments 
	 @param p probability of a success 
	 @return  quantil for prob probability 
	*/
	public final int InverseBinomial(double prob, int n, double p)
	{
		double acum = 0.0;
		for (int i = 0;i <= n;i++)
		{
			acum += Binomial(i, n, p);
			if (acum >= prob)
			{
				return i;
			}
		}
		return 0;
	}


	//endregion

	//region Statistics

	//region Mean

	/**  Mean of binomial distribution 
	 @param n number of independent yes/no experiments 
	 @param p probability of a success 
	 @return  mean value 
	*/
	public final double GetBinomialMean(int n, int p)
	{
		return n * p;
	}

	/**  Mean of hypergeometric distribution 
	 @param N number of independent yes/no experiments of the population 
	 @param n number of independent yes/no experiments of the sample 
	 @param d number of favorable elements in the sample 
	 @return  mean value 
	*/
	public final double GetHypergeometricMean(int N, int n, int d)
	{
		return (n * d) / N;
	}

	/**  Mean of poisson distribution 
	 @param lambda the expected number of occurrences that occur during the given interval 
	 @return  mean value 
	*/
	public final double GetPoissonMean(double lambda)
	{
		return lambda;
	}

	/**  Mean of negative binomial distribution distribution 
	 @param x number of experiments (fails) needed until n is reached 
	 @param n specified number of succeses 
	 @param p probability of a success 
	 @return  mean value 
	*/
	public final double GetNegativeBinomialMean(int x, int n, double p)
	{
		return x / p;
	}

	//endregion

	//region Standard Deviation

	/**  Standard deviation of binomial distribution 
	 @param n number of independent yes/no experiments 
	 @param p probability of a success 
	 @return  Standard deviation value 
	*/
	public final double GetBinomialStDev(int n, int p)
	{
		return Math.sqrt(n * p * (1 - p));
	}

	/**  Standard deviation of hypergeometric distribution 
	 @param N number of independent yes/no experiments of the population 
	 @param n number of independent yes/no experiments of the sample 
	 @param d number of favorable elements in the sample 
	 @return  Standard deviation value 
	*/
	public final double GetHypergeometricStDev(int N, int n, int d)
	{
		return (n * (d / N) * (1 - (d / N)) * (N - n)) / (N - 1);
	}

	/**  Standard deviation of poisson distribution 
	 @param lambda the expected number of occurrences that occur during the given interval 
	 @return  Standard deviation value 
	*/
	public final double GetPoissonStDev(double lambda)
	{
		return Math.sqrt(lambda);
	}

	/**  Standard deviation of negative binomial distribution distribution 
	 @param x number of experiments (fails) needed until n is reached 
	 @param n specified number of succeses 
	 @param p probability of a success 
	 @return  Standard deviation value 
	*/
	public final double GetNegativeBinomialStDev(int x, int n, double p)
	{
		return (x * (1 - p)) / Math.pow(p, 2);
	}

	//endregion

	//endregion

	//endregion

	//region Private Methods

	//region Distributions

	private double Binomial(int x, int n, double p)
	{
		return comb.Combinations(n, x) * Math.pow(p, x) * Math.pow((1 - p), (n - x));
	}


	private double Poisson(int x, double lambda)
	{
		return (Math.exp(-lambda) * Math.pow(lambda, x)) / comb.Factorial(x);
	}

	private double HurdlePoisson(int x, double alpha, double lambda)
	{
		if (x == 0)
		{
			return alpha;
		}
		return ((1 - alpha) * Math.pow(lambda, x)) / (comb.Factorial(x) * (Math.exp(lambda) - 1));
	}

	private double Hypergeometric(int x, int N, int n, int d)
	{
		return (comb.Combinations(d, x) * comb.Combinations(N - d, N - x)) / comb.Combinations(N, n);
	}

	private double NegativeBinomial(int x, int n, double p)
	{
		return comb.Combinations(x - 1, n - 1) * Math.pow(p, n) * Math.pow(1 - p, x - n);
	}

	//endregion

	// region Log Likelihood

	//n : number of periods - n0 : number of periods without event - t : total of events
	private double NegativeLogLikelihood(double alpha, double lambda, int n, int n0, int t)
	{
		return -n0 * Math.log(alpha) - (n - n0) * Math.log((1 - alpha) / (Math.exp(lambda) - 1)) - t * Math.log(lambda);
	}

	/**  Minimization based on likelihood 
	 @param alpha alpha parameter 
	 @param lambda lambda parameter (return parameter) 
	 @param maxLambda maximum for lambda 
	 @param n number of data 
	 @param n0 number of nonzero data 
	 @param t total of events 
	*/
	public final void MinimizeLikelihood(double alpha, DotNetHelpers.RefObject<Double> lambda, int maxLambda, int n, int n0, int t)
	{
		double minNll = Double.MAX_VALUE;
		double nll = 0;

		for (int l = 1; l < maxLambda; l++)
		{
			nll = NegativeLogLikelihood(alpha, (double)l, n, n0, t);
			if (nll < minNll)
			{
				lambda.argValue = (double)l;
				minNll = nll;
			}
		}

		double lam = lambda.argValue-1;
		double Lambda = lambda.argValue;
		while (lam <= lambda.argValue+1)
		{
			nll = NegativeLogLikelihood(alpha, lam, n, n0, t);
			if (nll < minNll)
			{
				Lambda = lam;
				minNll = nll;
			}
			lam = lam + 0.1;
		}
		lambda.argValue = Lambda;

		//region Obsolete

		/*
		   RandomGen rand = new RandomGen();
		   double a, l;

		   //explotarion: montecarlo
		   for (int i = 0; i < 100; i++) {
		       a = rand.Next(0, 99) / 100.0;
		       l = rand.Next(0, maxLambda);
		       TestParameters(a, l, ref alpha, ref lambda, n, n0, t, ref minNll);
		   }

		   //explotation: gradient
		   for (int i = 0; i < 100; i++) {
		       if (alpha > 0.01) { TestParameters(alpha - 0.01, lambda, ref alpha, ref lambda, n, n0, t, ref minNll); }
		       if (alpha < 0.99) { TestParameters(alpha + 0.01, lambda, ref alpha, ref lambda, n, n0, t, ref minNll); }
		       if (lambda > 1) { TestParameters(alpha, lambda - 1, ref alpha, ref lambda, n, n0, t, ref minNll); }
		       if (lambda < maxLambda) { TestParameters(alpha, lambda + 1, ref alpha, ref lambda, n, n0, t, ref minNll); }
		   }
		   */
		//endregion

	}

	private void TestParameters(double a, double l, DotNetHelpers.RefObject<Double> alpha, DotNetHelpers.RefObject<Double> lambda, int n, int n0, int t, DotNetHelpers.RefObject<Double> minNll)
	{
		double nll = NegativeLogLikelihood(a, l, n, n0, t);
		if (nll < minNll.argValue)
		{
			alpha.argValue = (double)a / 100.0;
			lambda.argValue = l;
			minNll.argValue = nll;
		}
	}

	//endregion

	//endregion

}