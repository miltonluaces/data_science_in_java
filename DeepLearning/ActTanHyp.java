package DeepLearning;

public class ActTanHyp extends ActFunction
{
	public final double Forward(double x)
	{
		return Math.tanh(x);
	}
	public final double Backward(double x)
	{
		double coshx = Math.cosh(x);
		double denom = (Math.cosh(2 * x) + 1);
		return 4 * coshx * coshx / (denom * denom);
	}
}