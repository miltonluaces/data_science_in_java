package SignalProcess;

import NumCalc.*;
import ClassicalStat.*;
import java.util.*;

//region Imports


//endregion


//region Comments

// Daubechies D4 wavelet transform (D4 denotes four coefficients)

//I have to confess up front that the comment here does not even come close to describing wavelet algorithms and the Daubechies D4
//algorithm in particular.  I don't think that it can be described in anything less than a journal article or perhaps a book.  I even have
//to apologize for the notation I use to describe the algorithm, which is barely adequate.  But explaining the correct notation would take
//a fair amount of space as well.  This comment really represents some notes that I wrote up as I implemented the code.  If you are
//unfamiliar with wavelets I suggest that you look at the bearcave.com web pages and at the wavelet literature.  I have yet to see a really
//good reference on wavelets for the software developer.  The best book I can recommend is <i>Ripples in Mathematics</i> by Jensen and Cour-Harbo.
//All wavelet algorithms have two components, a wavelet function and a  scaling function.  These are sometime also referred to as high pass
//and low pass filters respectively.
//The wavelet function is passed two or more samples and calculates a wavelet coefficient.  In the case of the Haar wavelet this is 
//The scaling function produces a smoother version of the original data.  In the case of the Haar wavelet algorithm
//this is an average of two adjacent elements.
//The Daubechies D4 wavelet algorithm also has a wavelet and a scaling function.  The coefficients for the scaling function are denoted as h<sub>i</sub> and the
//wavelet coefficients are g<sub>i</sub>.
//Mathematicians like to talk about wavelets in terms of   a wavelet algorithm applied to an infinite data set.
//In this case one step of the forward transform can be expressed  as the infinite matrix of wavelet coefficients
//represented below multiplied by the infinite signal vector.

//   a<sub>i</sub> = ...h0,h1,h2,h3, 0, 0, 0, 0, 0, 0, 0, ...   s<sub>i</sub>
//   c<sub>i</sub> = ...g0,g1,g2,g3, 0, 0, 0, 0, 0, 0, 0, ...   s<sub>i+1</sub>
//  a<sub>i+1</sub> = ...0, 0, h0,h1,h2,h3, 0, 0, 0, 0, 0, ...   s<sub>i+2</sub>
//  c<sub>i+1</sub> = ...0, 0, g0,g1,g2,g3, 0, 0, 0, 0, 0, ...   s<sub>i+3</sub>
//  a<sub>i+2</sub> = ...0, 0, 0, 0, h0,h1,h2,h3, 0, 0, 0, ...   s<sub>i+4</sub>
//  c<sub>i+2</sub> = ...0, 0, 0, 0, g0,g1,g2,g3, 0, 0, 0, ...   s<sub>i+5</sub>
//  a<sub>i+3</sub> = ...0, 0, 0, 0, 0, 0, h0,h1,h2,h3, 0, ...   s<sub>i+6</sub>
//  c<sub>i+3</sub> = ...0, 0, 0, 0, 0, 0, g0,g1,g2,g3, 0, ...   s<sub>i+7</sub>

 //       The dot product (inner product) of the infinite vector and
 // a row of the matrix produces either a smoother version of the
 // signal (a<sub>i</sub>) or a wavelet coefficient (c<sub>i</sub>).
 // </p>
 // <p>
 // In an ordered wavelet transform, the smoothed (a<sub>i</sub>) are
 // stored in the first half of an <i>n</i> element array region.  The
 // wavelet coefficients (c<sub>i</sub>) are stored in the second half
 // the <i>n</i> element region.  The algorithm is recursive.  The
 // smoothed values become the input to the next step.
 // </p>
 // <p>
 // The transpose of the forward transform matrix above is used
 // to calculate an inverse transform step.  Here the dot product is
 // formed from the result of the forward transform and an inverse
 // transform matrix row.
 // </p>
 // <pre>
 //     s<sub>i</sub> = ...h2,g2,h0,g0, 0, 0, 0, 0, 0, 0, 0, ...  a<sub>i</sub>
 //   s<sub>i+1</sub> = ...h3,g3,h1,g1, 0, 0, 0, 0, 0, 0, 0, ...  c<sub>i</sub>
 //   s<sub>i+2</sub> = ...0, 0, h2,g2,h0,g0, 0, 0, 0, 0, 0, ...  a<sub>i+1</sub>
 //   s<sub>i+3</sub> = ...0, 0, h3,g3,h1,g1, 0, 0, 0, 0, 0, ...  c<sub>i+1</sub>
 //   s<sub>i+4</sub> = ...0, 0, 0, 0, h2,g2,h0,g0, 0, 0, 0, ...  a<sub>i+2</sub>
 //   s<sub>i+5</sub> = ...0, 0, 0, 0, h3,g3,h1,g1, 0, 0, 0, ...  c<sub>i+2</sub>
 //   s<sub>i+6</sub> = ...0, 0, 0, 0, 0, 0, h2,g2,h0,g0, 0, ...  a<sub>i+3</sub>
 //   s<sub>i+7</sub> = ...0, 0, 0, 0, 0, 0, h3,g3,h1,g1, 0, ...  c<sub>i+3</sub>
 // </pre>

 // <p>
 // Using a standard dot product is grossly inefficient since most
 // of the operands are zero.  In practice the wavelet coefficient 
 // values are moved along the signal vector and a four element 
 // dot product is calculated.  Expressed in terms of arrays, for
 // the forward transform this would be:
 // </p>
 // <pre>
 // a<sub>i</sub> = s[i]*h0 + s[i+1]*h1 + s[i+2]*h2 + s[i+3]*h3
 // c<sub>i</sub> = s[i]*g0 + s[i+1]*g1 + s[i+2]*g2 + s[i+3]*g3
 // </pre>
 // <p>
 // This works fine if we have an infinite data set, since we don't
 // have to worry about shifting the coefficients "off the end" of
 // the signal.
 // </p>
 // <p>
 // I sometimes joke that I left my infinite data set in my other bear
 // suit.  The only problem with the algorithm described so far is that
 // we don't have an infinite signal.  The signal is finite.  In fact
 // not only must the signal be finite, but it must have a power of two
 // number of elements.
 // </p>
 // <p>
 // If i=N-1, the i+2 and i+3 elements will be beyond the end of 
 // the array.  There are a number of methods for handling the 
 // wavelet edge problem.  This version of the algorithm acts 
 // like the data is periodic, where the data at the start of 
 // the signal wraps around to the end.
 // </p>
 // <p>
 // This algorithm uses a temporary array.  A Lifting Scheme version of
 // the Daubechies D4 algorithm does not require a temporary.  The
 // matrix discussion above is based on material from <i>Ripples in
 // Mathematics</i>, by Jensen and Cour-Harbo.  Any error are mine.
 // </p>

 // <p>
 // <b>Author</b>: Ian Kaplan<br>
 // <b>Use</b>: You may use this software for any purpose as long
 // as I cannot be held liable for the result.  Please credit me
 // with authorship if use use this source code.
 // </p>
//endregion

/**  Algorithmic repository class for Daubechies Wavelet calculation 
*/
public class Daubechies
{

	//region Fields

	//scaling coeffs
	private double h0;
	private double h1;
	private double h2;
	private double h3;

	//wavelet coeffs
	private double g0;
	private double g1;
	private double g2;
	private double g3;

	//endregion

	//region Constructor

	/**  Constructor 
	*/
	public Daubechies()
	{

		double sqrt2 = Math.sqrt(2);
		double sqrt3 = Math.sqrt(3);

		h0 = (1 + sqrt3) / (4 * sqrt2);
		h1 = (3 + sqrt3) / (4 * sqrt2);
		h2 = (3 - sqrt3) / (4 * sqrt2);
		h3 = (1 - sqrt3) / (4 * sqrt2);

		g0 = h3;
		g1 = -h2;
		g2 = h1;
		g3 = -h0;
	}

	//endregion

	//region Public Methods

	   //Forward Daubechies D4 transform
	/**  Daubechies tranformed function 
	 @param serie time series 
	 @return  transformed time series 
	*/
	public final ArrayList<Double> daubTrans(ArrayList<Double> serie)
	{
		ArrayList<Double> trans = new ArrayList<Double>(serie);
		for (int n = serie.size();n >= 4;n = n / 2)
		{
			transform(trans, n);
		}
		return trans;
	}


	/**  Inverse Daubechies D4 transform 
	 @param coeff list of coefficents 
	 @return  list of calculated inverse values 
	*/
	public final ArrayList<Double> invDaubTrans(ArrayList<Double> coeff)
	{
		ArrayList<Double> coeffs = new ArrayList<Double>(coeff);
		for (int n = 4;n <= coeff.size();n = n * 2)
		{
			invTransform(coeffs, n);
		}
		return coeffs;
	}

	//endregion

	//region Private Methods

		//Forward wavelet transform. 
		private void transform(ArrayList<Double> serie, int n)
		{
			if (n < 4)
			{
				return;
			}
			double[] trans = new double[n];
			int i;

			//funciones D4 scaling y wavelet
			for (i = 0;i <= (n - 4) / 2;i++)
			{
				trans[i] = h0 * serie.get(2 * i) + h1 * serie.get(2 * i + 1) + h2 * serie.get(2 * i + 2) + h3 * serie.get(2 * i + 3); //scaling
				trans[i + n / 2] = g0 * serie.get(2 * i) + g1 * serie.get(2 * i + 1) + g2 * serie.get(2 * i + 2) + g3 * serie.get(2 * i + 3); //wavelet
			}
			//iteración final: como no hay n, n+1 se toma 0 y 1 para completar
			trans[i] = h0 * serie.get(n - 2) + h1 * serie.get(n - 1) + h2 * serie.get(0) + h3 * serie.get(1);
			trans[i + n / 2] = g0 * serie.get(n - 2) + g1 * serie.get(n - 1) + g2 * serie.get(0) + g3 * serie.get(1);

			//se copian valores en la serie original
			for (i = 0;i < n;i++)
			{
				serie.set(i, trans[i]);
			}

		}


		private void invTransform(ArrayList<Double> coeff, int n)
		{
			if (n >= 4)
			{

				double[] invTrans = new double[n];

				//iteracion inicial: como no hay -1 y -2, se toma ult y penultimo
				invTrans[0] = h2 * coeff.get(n / 2 - 1) + g2 * coeff.get(n - 1) + h0 * coeff.get(0) + g0 * coeff.get(n / 2);
				invTrans[1] = h3 * coeff.get(n / 2 - 1) + g3 * coeff.get(n - 1) + h1 * coeff.get(0) + g1 * coeff.get(n / 2);

				for (int i = 0;i < n / 2 - 1;i++)
				{
					invTrans[2 * (i + 1)] = h2 * coeff.get(i) + g2 * coeff.get(i + n / 2) + h0 * coeff.get(i + 1) + g0 * coeff.get(i + n / 2 + 1); //scaling
					invTrans[2 * (i + 1) + 1] = h3 * coeff.get(i) + g3 * coeff.get(i + n / 2) + h1 * coeff.get(i + 1) + g1 * coeff.get(i + n / 2 + 1); //wavelet
				}

				//se copian valores en la serie original
				for (int i = 0;i < n;i++)
				{
					coeff.set(i, invTrans[i]);
				}
			}
		}


//endregion

}