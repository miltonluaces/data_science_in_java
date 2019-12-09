package DeepLearning;

//region Imports


//endregion

public abstract class ActFunction
{

	//region Fields

	private static long uid = 1;
	private long id;

	//endregion

	//region Constructor

	public ActFunction()
	{
		id = uid + 1;
	}

	//endregion

	//region Properties
	public final long getId()
	{
		return id;
	}

	//endregion

	//region Public Methods
	public double Forward(double x)
	{
		return 0;
	}
	public double Backward(double x)
	{
		return 0;
	}

	//endregion
}