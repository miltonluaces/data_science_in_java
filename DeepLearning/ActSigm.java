package DeepLearning;

public class ActSigm extends ActFunction
{
	@Override
	public double Forward(double x)
	{
		return 1 / (1 + Math.exp(-x));
	}
	@Override
	public double Backward(double x)
	{
		double act = Forward(x);
		return act * (1 - act);
	}
}