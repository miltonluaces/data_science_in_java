package NumCalc;

import java.util.*;

/**  Repository of derivative calculation methods 
*/
public class Derivatives
{

	//region Constructor 

	/**  Constructor  
	*/
	public Derivatives()
	{

	}

	//endregion 

	//region Public Methods

	/**  Get the derivative of any order, of an array of data (Y) with the corresponding independent variable (X) 
	 @param X independent variable 
	 @param Y data (dependent variable) 
	 @param order derivative order 
	 @return  list of calculated data 
	*/
	public final ArrayList<Double> GetDerivative(ArrayList<Double> X, ArrayList<Double> Y, int order)
	{
		if (X.size() != Y.size())
		{
			throw new RuntimeException("X_and_Y_must_have_the_same_size");
		}

		ArrayList<Double> D = Y;
		for (int i = 0;i < order;i++)
		{
			D = GetDerivative(X, D);
		}
		return D;
	}

	/**  Get the derivative of any order, of an array of data (Y) 
	 @param Y data 
	 @param order derivative order 
	 @return  list of calculated data 
	*/
	public final ArrayList<Double> GetDerivative(ArrayList<Double> Y, int order)
	{

		ArrayList<Double> D = Y;
		for (int i = 0;i < order;i++)
		{
			D = GetDerivative(D);
		}
		return D;
	}

	//endregion 

	//region Private Methods

	private ArrayList<Double> GetDerivative(ArrayList<Double> X, ArrayList<Double> Y)
	{
		if (X.size() != Y.size())
		{
			throw new RuntimeException("X_and_Y_must_have_the_same_size");
		}

		ArrayList<Double> deriv = new ArrayList<Double>();
		double diffQuotient = 0.0;
		for (int i = 0;i < X.size() - 1;i++)
		{
			diffQuotient = (Y.get(i + 1) - Y.get(i)) / (X.get(i + 1) - X.get(i));
			deriv.add(diffQuotient);
		}
		deriv.add(diffQuotient);
		return deriv;
	}

	private ArrayList<Double> GetDerivative(ArrayList<Double> Y)
	{

		ArrayList<Double> deriv = new ArrayList<Double>();
		double diffQuotient = 0.0;
		for (int i = 0;i < Y.size() - 1;i++)
		{
			diffQuotient = (Y.get(i + 1) - Y.get(i));
			deriv.add(diffQuotient);
		}
		deriv.add(diffQuotient);
		return deriv;
	}

	//endregion

}