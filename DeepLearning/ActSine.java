package DeepLearning;

public class ActSine extends ActFunction
{
	public final double Forward(double x)
	{
		return Math.sin(x);
	}
	public final double Backward(double x)
	{
		return Math.cos(x);
	}
}