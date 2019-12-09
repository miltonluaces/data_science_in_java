package NumCalc;

//region 


//endregion


public class MemoryCalc
{

	private double x1;
	private double x2;
	private double y1;
	private double y2;

	private double a;
	private double b;

	public final void CalcAB(double x1, double x2, double y1, double y2)
	{
		a = (y2 - y1) / (x2 - x1);
		b = y1 - (x1 * (y2 - y1) / (x2 - x1));
	}

	public final double LineEquation(double x)
	{
		return a * x + b;
	}

	public final double InverseLineEquation(double y)
	{
		return (y - b) / a;
	}
}