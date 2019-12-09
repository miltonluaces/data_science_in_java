package NumCalc;

//region Imports


//endregion


public class KalmanFilter
{

	//region Fields

	private Matrix x0;
	private Matrix p0;

	private Matrix f;
	private Matrix b;
	private Matrix u;
	private Matrix q;
	private Matrix h;
	private Matrix r;
	private Matrix state;
	private Matrix cov;

	//endregion

	//region Constructor

	 public KalmanFilter(Matrix f, Matrix b, Matrix u, Matrix q, Matrix h, Matrix r)
	 {
		this.f = f;
		this.b = b;
		this.u = u;
		this.q = q;
		this.h = h;
		this.r = r;
	 }

	//endregion

	//region Properties

	public final Matrix getF()
	{
		return f;
	}
	public final Matrix getB()
	{
		return b;
	}
	public final Matrix getU()
	{
		return u;
	}
	public final Matrix getQ()
	{
		return q;
	}
	public final Matrix getH()
	{
		return h;
	}
	public final Matrix getR()
	{
		return r;
	}
	public final Matrix getState()
	{
		return state;
	}
	public final Matrix getCov()
	{
		return cov;
	}
	public final Matrix getX0()
	{
		return x0;
	}
	public final Matrix getP0()
	{
		return p0;
	}

	//endregion

	//region Public Mehtod

	/*
	public void Predict() {
	    x0 = f.Multiply(state) + (b.Multiply(u));
	    p0 = f.Multiply(cov).Multiply(f.Transpose()) + q;
	}

	public void Correct(Matrix z) {
	    Matrix s = h.Multiply(p0) * h.Transpose() + r;
	    Matrix k = p0.Multiply(h).Transpose() * s.Inverse();
	    state = x0 + (k * (z - (h.Multiply(x0)));
	    Matrix I = Matrix.Identity(p0.Rows, p0.Columns);
	    cov = (I - k.Multiply(h)) * p0;
	}
	*/

	//endregion
}