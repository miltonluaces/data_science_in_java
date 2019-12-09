package SignalProcess;

//region Imports


//endregion


/**  Data class for complex numbers 
*/
public class Complex
{

	//region Fields

	private double real;
	private double imag;

	//endregion

	//region Constructor

	/**  Constructor 
	 @param real real component 
	 @param imag imaginary component 
	*/
	public Complex(double real, double imag)
	{
		this.real = real;
		this.imag = imag;
	}

	//endregion

	//region Properties

	/**  real part 
	*/
	public final double getReal()
	{
		return real;
	}

	/**  imaginary part 
	*/
	public final double getImag()
	{
		return imag;
	}

	/**  abs/modulus/magnitud 
	*/

	public final double getAbs()
	{
		return Math.sqrt(real * real + imag * imag);
	}

	/**  angle/phase/argument  between -pi and pi 
	*/
	public final double getPhase()
	{
		return Math.atan2(imag, real);
	}

	//endregion

	//region Operations

	/**  Plus: return a new Complex object whose value is (this + b) 
	 @param b the complex number for the operation 
	 @return  complex number result 
	*/
	public final Complex Plus(Complex b)
	{
		Complex a = this;
		double real = a.real + b.real;
		double imag = a.imag + b.imag;
		return new Complex(real, imag);
	}

	/**  Minus:  return a new Complex object whose value is (this - b) 
	 @param b the complex number for the operation 
	 @return  complex number result 
	*/
	public final Complex Minus(Complex b)
	{
		Complex a = this;
		double real = a.real - b.real;
		double imag = a.imag - b.imag;
		return new Complex(real, imag);
	}

	/**  Times:  return a new Complex object whose value is (this * b) 
	 @param b the complex number for the operation 
	 @return  complex number result 
	*/
	public final Complex Times(Complex b)
	{
		Complex a = this;
		double real = a.real * b.real - a.imag * b.imag;
		double imag = a.real * b.imag + a.imag * b.real;
		return new Complex(real, imag);
	}

	/**  Scalar multiplication: returns a new object whose value is (this * alpha) 
	 @param alpha factor 
	 @return  complex number result 
	*/
	public final Complex Times(double alpha)
	{
		return new Complex(alpha * real, alpha * imag);
	}

	/**  Conjugate: returns a new Complex object whose value is the conjugate of this 
	 @return  complex number result 
	*/
	public final Complex Conjugate()
	{
		return new Complex(real, -imag);
	}

	/**  Reciprocal: returns a new Complex object whose value is the reciprocal of this 
	 @return  complex number result 
	*/
	public final Complex Reciprocal()
	{
		double scale = real * real + imag * imag;
		return new Complex(real / scale, -imag / scale);
	}

	/**  Divides: returns a / b 
	 @param b divisor 
	 @return  complex number result 
	*/
	public final Complex Divides(Complex b)
	{
		Complex a = this;
		return a.Times(b.Reciprocal());
	}

	/**  Exp: returns a new Complex object whose value is the complex exponential of this 
	 @return  complex number result 
	*/
	public final Complex Exp()
	{
		return new Complex(Math.exp(real) * Math.cos(imag), Math.exp(real) * Math.sin(imag));
	}

	/**  Sin: returns a new Complex object whose value is the complex sine of this 
	 @return  complex number result 
	*/
	public final Complex Sin()
	{
		return new Complex(Math.sin(real) * Math.cosh(imag), Math.cos(real) * Math.sinh(imag));
	}

	/**  Cos: returns a new Complex object whose value is the complex cosine of this 
	 @return  complex number result 
	*/
	public final Complex Cos()
	{
		return new Complex(Math.cos(real) * Math.cosh(imag), -Math.sin(real) * Math.sinh(imag));
	}

	/**  Tan: returns a new Complex object whose value is the complex tangent of this 
	 @return  complex number result 
	*/
	public final Complex Tan()
	{
		return Sin().Divides(Cos());
	}

	//Static Plus
	/**  Sum of two numbers, static 
	 @param a first number 
	 @param b second number 
	 @return  complex number result 
	*/
	public static Complex Plus(Complex a, Complex b)
	{
		double real = a.real + b.real;
		double imag = a.imag + b.imag;
		Complex sum = new Complex(real, imag);
		return sum;
	}

	//endregion

	//region Override To String

	/**  ToString override  
	 @return  the formatted string 
	*/
	@Override
	public String toString()
	{
		if (imag == 0)
		{
			return real + "";
		}
		if (real == 0)
		{
			return imag + " i";
		}
		if (imag < 0)
		{
			return real + " - " + (-imag) + " i";
		}
		return real + " + " + imag + " i";
	}

	//endregion

}