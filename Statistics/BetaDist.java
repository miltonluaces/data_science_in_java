package ClassicalStat;

//region Imports


//endregion


/**  Beta Distribution 
*/
public class BetaDist
{

	//region Fields

	private Statistics stat;
	private double a;
	private double b;


	//endregion

	//region Constructor

	/**  Constructor 
	*/
	public BetaDist() {
		stat = new Statistics();
	}

	//endregion

	//region Properties

	/**  a parameter 
	*/
	public final double getA() {
		return a;
	}
	public final void setA(double value) {
		a = value;
	}

	/**  b parameter 
	*/
	public final double getB() 	{
		return b;
	}
	public final void setB(double value) {
		b = value;
	}


	//endregion

	//region Public Methods

	/**  Probability function of a parameter theta 
	 @param theta the value of the parameter 
	 @param acum if accumulated or not 
	 @return  the probability 
	*/
	public final double Probability(double theta, boolean acum)
	{
		if (acum)
		{
			return ProbabilityAcum(theta);
		}
		else
		{
			return Probability(theta);
		}
	}

	private double Probability(double theta)
	{
		if (theta < 0 || theta > 1)
		{
			return 0;
		}
		double num = Math.pow(theta, (a - 1)) * Math.pow((1 - theta),(b - 1));
		double den = stat.B(a, b);
		return num / den;
	}

	private double ProbabilityAcum(double theta)
	{
		return stat.BetaInc(a, b, theta) / stat.B(a, b);
	}

	/**  Quantile 
	 @param p the probability 
	 @return  the quantile 
	*/
	public final double Quantile(double p)
	{
		return MonotoneBisection((double x) -> ProbabilityAcum(x), true, 0.0, 1.0, p, 0.01);
	}

	/**  Probability function delegate 
	 @param x value 
	 @return  probability 
	*/
	@FunctionalInterface
	public interface function
	{
		double invoke(double x);
	}

	/**  Monotone bisection for beta distribtuion 
	 @param f Probability function 
	 @param ascendant if it is ascendant or not 
	 @param min minimum value 
	 @param max maximum value 
	 @param val value to obtain the inverse 
	 @param eps epsilon value 
	 @return  the inverse value 
	*/
	public final double MonotoneBisection(function f, boolean ascendant, double min, double max, double val, double eps)
	{
		int it = 100;
		int maxIt = 100;
		DotNetHelpers.RefObject<Integer> tempRef_it = new DotNetHelpers.RefObject<Integer>(it);
		double tempVar = MonotoneBisectionRec(f, min, max, val, ascendant, eps, tempRef_it, maxIt);
		it = tempRef_it.argValue;
		return tempVar;
	}

	private double MonotoneBisectionRec(function f, double min, double max, double val, boolean ascendant, double eps, DotNetHelpers.RefObject<Integer> it, int maxIt)
	{
		it.argValue++;
		double inter = min + (max - min) / 2.0;
		if (it.argValue > maxIt || inter == min || inter == max)
		{
			//Console.WriteLine("Error.More than 30 iterations.");
			return inter;
		}
		double fInter = f.invoke(inter);
		if (Math.abs(val - fInter) < eps)
		{
			return inter;
		}
		if (ascendant)
		{
			if (val < fInter)
			{
				return MonotoneBisectionRec(f, min, inter, val, ascendant, eps, it, maxIt);
			}
			if (val > fInter)
			{
				return MonotoneBisectionRec(f, inter, max, val, ascendant, eps, it, maxIt);
			}
		}
		else
		{
			if (val > fInter)
			{
				return MonotoneBisectionRec(f, min, inter, val, ascendant, eps, it, maxIt);
			}
			if (val < fInter)
			{
				return MonotoneBisectionRec(f, inter, max, val, ascendant, eps, it, maxIt);
			}
		}
		//Console.WriteLine("Error inter = " + inter);
		return inter;
	}

	/**  Expectation 
	 @return  the mean 
	*/
	public final double Mean()
	{
		return a / (a + b);
	}

	/**  Variance 
	 @return  the variance 
	*/
	public final double Var()
	{
		return (a * b) / (Math.pow(a + b, 2) * (a + b + 1));
	}

	/**  Set B parameter having a parameter and the mean (solving for b) for elicitation 
	 @param a a parameter 
	 @param m expectation (mean) 
	*/
	public final void SetB(int a, double m)
	{
		b = a * (1 - m) / m;
	}
}

	//endregion

