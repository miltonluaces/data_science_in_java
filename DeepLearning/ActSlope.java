package DeepLearning;

public class ActSlope extends ActFunction
{

	private double slope = 0;
	public ActSlope(double slope)
	{
		super();
	this.slope = slope;
	}

	public final double Forward(double x)
	{
		return (x >= 0) ? x : x * slope;
	}
	public final double Backward(double x)
	{
		return (x >= 0) ? 1 : slope;
	}
}