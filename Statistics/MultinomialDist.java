package ClassicalStat;

import java.util.*;

//region Imports


//endregion


/**  Multinomial distribution class 
*/
public class MultinomialDist
{

	//region Fields

	private ArrayList<Double> P;
	private Combinatory comb;

	//endregion

	//region Constructor

	/**  Constructor 
	*/
	public MultinomialDist()
	{
		P = new ArrayList<Double>();
		comb = new Combinatory();
	}

	//endregion

	//region Public Methods

	/**  Set proportions 
	 @param P vector of proportions (sum 1) 
	*/
	public final void SetP(ArrayList<Double> P)
	{
		double sum = 0;
		for (double p : P)
		{
			sum += p;
		}
		if (sum != 1)
		{
			throw new RuntimeException("Error. Probabilities should sum to 1");
		}
		this.P = P;
	}

	/**  Probability of a specific X value 
	 @param X vector X value 
	 @return  the probability 
	*/
	public final double Probability(List<Integer> X)
	{
		int n = 0;
		for (int x : X)
		{
			n += x;
		}
		double factNum = comb.Factorial(n);
		double prodPExp = 1;
		double prodXFact = 1;
		for (int i = 0;i < P.size();i++)
		{
			prodPExp *= Math.pow(P.get(i), X.get(i));
			prodXFact *= comb.Factorial(X.get(i));
		}
		return (factNum / prodXFact) * prodPExp;
	}

	//endregion

}