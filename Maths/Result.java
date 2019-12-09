package NumCalc;

import ClassicalStat.*;
import java.util.*;

//region Inner class Result

/**  Result of a Regression calculation 
*/
public class Result
{

	//region Fields

	/**  first coefficent 
	*/
	public double c0;

	/**  second coefficent 
	*/
	public double c1;

	/**  following first coefficent 
	*/
	public double sig0;

	/**  following second coefficent 
	*/
	public double sig1;

	/**  Squared chi coefficent 
	*/
	public double chi2;

	/**  Determination coefficent 
	*/
	public double r2;

	/**  list of weights 
	*/
	public ArrayList<Double> W;

	private int n;

	//endregion

	//region Constructors

	/**  Constructor 
	 @param n number of data 
	*/
	public Result(int n)
	{
		this.n = n;
		this.c0 = 0;
		this.c1 = 0;
		this.sig0 = 0;
		this.sig1 = 0;
		this.chi2 = 0;
		this.r2 = 0;
		this.W = new ArrayList<Double>();
		for (int i = 0;i < n;i++)
		{
			W.add(1.0);
		}
	}

	/**  Constructor 
	 @param n number of data 
	 @param W list of weights 
	*/
	public Result(int n, ArrayList<Double> W)
	{
		this.n = n;
		this.c0 = 0;
		this.c1 = 0;
		this.sig0 = 0;
		this.sig1 = 0;
		this.chi2 = 0;
		this.r2 = 0;
		this.W = W;
	}

	//endregion

	//region Public Methods

	/**  Resets all weights to 1.0 (to avoid weighting effect) 
	*/
	public final void ResetW()
	{
		for (int i = 0;i < n;i++)
		{
			W.set(i, 1.0);
		}
	}

	//endregion
}
//endregion
