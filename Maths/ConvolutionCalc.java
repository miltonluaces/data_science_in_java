package NumCalc;

//region Imports

import ClassicalStat.*;
import java.util.*;

//endregion


/**  Interface for all convolution calculators 
*/
public interface ConvolutionCalc
{

	/**  Load data for calculation 
	 @param data data from time series 
	 @param n number of convolutions 
	*/
	void LoadData(ArrayList<Double> data, double n);

	/**  Load histogram for calculation 
	 @param hist histogram for calculation 
	 @param n number of convolutions 
	*/
	void LoadHistogram(Histogram hist, double n);

	/**  Acumulated probability for a certain value 
	 @param x the value 
	 @return  the probability 
	*/
	double ProbabilityAcum(double x);

	/**  Quantile for a certain probability 
	 @param p the probability 
	 @return  the value 
	*/
	double Quantile(double p);

	/**  if calculation is valid 
	 @return  true if it is valid, false if not 
	*/
	boolean Valid();


}