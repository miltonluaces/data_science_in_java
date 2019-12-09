package SignalProcess;

import java.util.*;

/**  Fast Fourier Transform and Inverse Fast Fourier Transform (N = power of 2.) 
*/
public class FFT extends SignalProc
{

	//region Constructor

	/**  Constructor 
	*/
	public FFT()
	{
	}

	//endregion

	//region Real List Methods

	/**  Fast fourier transform 
	 @param X list of original data (independent vector variable) 
	 @return  list of transformed data 
	*/
	public final ArrayList<Double> FftReg(ArrayList<Double> X)
	{
		Complex[] XComp = new Complex[X.size()];
		for (int i = 0;i < X.size();i++)
		{
			XComp[i] = new Complex(X.get(i), 0);
		}
		Complex[] YComp = Fft(XComp);
		XComp = Ifft(YComp);
		ArrayList<Double> Y = new ArrayList<Double>();
		for (Complex xComp : XComp)
		{
			Y.add(xComp.Real);
		}
		return Y;
	}

	 //endregion

	//region Complex Array Methods

	/**  Fast Fourier Transform : Cooley-Tukey algorithm
	 @param X list of complex values 
	 @return  array of complex transformed data 
	*/
	public final Complex[] Fft(Complex[] X)
	{
		int N = X.length;

		if (N == 1)
		{
			return new Complex[] {X[0]};
		}
		if (N % 2 != 0)
		{
			throw new RuntimeException(Strings.N_is_not_a_power_of_2);
		}

		//generate fft of even and odd integers
		Complex[] even = new Complex[N / 2];
		for (int k = 0;k < N / 2;k++)
		{
			even[k] = X[2 * k];
		}
		Complex[] fftEven = Fft(even);
		Complex[] odd = new Complex[N / 2];
		for (int k = 0;k < N / 2;k++)
		{
			odd[k] = X[2 * k + 1];
		}
		Complex[] fftOdd = Fft(odd);

		//combine terms
		Complex[] Y = new Complex[N];
		for (int k = 0;k < N / 2;k++)
		{
			double kth = -2 * k * Math.PI / N;
			Complex wk = new Complex(Math.cos(kth), Math.sin(kth));
			Y[k] = fftEven[k].Plus(wk.Times(fftOdd[k]));
			Y[k + N / 2] = fftEven[k].Minus(wk.Times(fftOdd[k]));
		}
		return Y;
	}


	/**  Inverse Fast Fourier Transform 
	 @param Y transformed function 
	 @return  Complex transformed array 
	*/
	public final Complex[] Ifft(Complex[] Y)
	{
		int N = Y.length;
		Complex[] X = new Complex[N];

		//take conjugate
		for (int i = 0;i < N;i++)
		{
			X[i] = Y[i].Conjugate();
		}

		//compute forward FFT
		X = Fft(X);

		//take conjugate again
		for (int i = 0;i < N;i++)
		{
			X[i] = X[i].Conjugate();
		}

		//divide by N
		for (int i = 0;i < N;i++)
		{
			X[i] = X[i].Times(1.0 / N);
		}

		return X;
	}

	/**  Circular convolution of x and y 
	 @param Χ list of complex independent variable vector data 
	 @param Υ  list of complex dependent variable vector data 
	 @return  Complex array of the convolution 
	*/
	public final Complex[] CircularConvolve(Complex[] Χ, Complex[] Υ)
	{

		//should probably pad x and y with 0s so that they have same length and are powers of 2
		if (Χ.length != Υ.length)
		{
			throw new RuntimeException(Strings.Dimensions_do_not_agree);
		}

		int N = Χ.length;

		//compute FFT of each sequence
		Complex[] A = Fft(Χ);
		Complex[] B = Fft(Υ);

		//point-wise multiply
		Complex[] C = new Complex[N];
		for (int i = 0;i < N;i++)
		{
			C[i] = A[i].Times(B[i]);
		}

		//compute inverse FFT
		return Ifft(C);
	}


	//Linear convolution of x and y
	/**  Linear convolution 
	 @param Χ list of complex independent variable vector data
	 @param Y list of complex dependent variable vector data
	 @return  Complex array of the linear convolution 
	*/
	public final Complex[] LinearConvolve(Complex[] Χ, Complex[] Y)
	{
		Complex zero = new Complex(0, 0);

		Complex[] A = new Complex[2 * Χ.length];
		for (int i = 0;i < Χ.length;i++)
		{
			A[i] = Χ[i];
		}
		for (int i = Χ.length;i < 2 * Χ.length;i++)
		{
			A[i] = zero;
		}

		Complex[] B = new Complex[2 * Y.length];
		for (int i = 0;i < Y.length;i++)
		{
			B[i] = Y[i];
		}
		for (int i = Y.length;i < 2 * Y.length;i++)
		{
			B[i] = zero;
		}

		return CircularConvolve(A, B);
	}

	/**  Circular convolution of x and y 
	 @param Χ list of complex independent variable vector data 
	 @return  Complex array of the convolution 
	*/
	public final Complex[] CircularConvolve(Complex[] Χ)
	{

		int N = Χ.length;

		//compute FFT of the sequence
		Complex[] A = Fft(Χ);

		//point-wise multiply
		Complex[] C = new Complex[N];
		for (int i = 0;i < N;i++)
		{
			C[i] = A[i].Times(A[i]);
		}

		//compute inverse FFT
		return Ifft(C);
	}

	//Linear convolution of x and y
	/**  Linear convolution 
	 @param Χ list of complex independent variable vector data
	 @return  Complex array of the linear convolution 
	*/
	public final Complex[] LinearConvolve(Complex[] Χ)
	{
		Complex zero = new Complex(0, 0);

		Complex[] A = new Complex[2 * Χ.length];
		for (int i = 0;i < Χ.length;i++)
		{
			A[i] = Χ[i];
		}
		for (int i = Χ.length;i < 2 * Χ.length;i++)
		{
			A[i] = zero;
		}

		return CircularConvolve(A);
	}

	/**  Multiple fft convolution 
	 @param Data data complex array 
	 @param n number of convolution 
	 @return  Complex array of result 
	*/
	public final Complex[] LinearConvolve(Complex[] Data, int n)
	{
		Complex zero = new Complex(0, 0);
		Complex two = new Complex(2, 0);
		Complex[] A = null;
		int length = Data.length;
		Complex[] C = null;
		Complex[] CAnt = new Complex[Data.length * 2];
		for (int i = 0;i < Data.length;i++)
		{
			CAnt[i] = Data[i];
		}
		for (int i = Data.length;i < Data.length * 2;i++)
		{
			CAnt[i] = zero;
		}
		for (int i = 0;i < n;i++)
		{
			length *= 2;
			A = new Complex[length];
			for (int j = 0;j < Data.length;j++)
			{
				A[j] = Data[j];
			}
			for (int j = Data.length;j < length;j++)
			{
				A[j] = zero;
			}
			if (C != null)
			{
				CAnt = C;
			}
			C = new Complex[length];
			for (int j = 0;j < CAnt.length;j++)
			{
				C[j] = CAnt[j];
			}
			for (int j = CAnt.length;j < length;j++)
			{
				C[j] = zero;
			}
			C = CircularConvolve(C, A);

		}
		return C;
	}

	//endregion
}