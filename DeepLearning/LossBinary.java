package DeepLearning;

//region Imports


//endregion

public class LossBinary implements ILossFunction
{

	//region ILossFunction Implementation

	public final double Measure(Neuron act, Neuron exp)
	{
		if (act.W.getLength() != exp.W.getLength())
		{
			throw new RuntimeException("mismatch");
		}

		for (int i = 0; i < exp.W.getLength(); i++)
		{
			if (exp.W[i] >= 0.5 && act.W[i] < 0.5)
			{
				return 1;
			}
			if (exp.W[i] < 0.5 && act.W[i] >= 0.5)
			{
				return 1;
			}
		}
		return 0;
	}

	public final void Backward(Neuron act, Neuron exp)
	{
		throw new UnsupportedOperationException("not implemented");
	}

	//endregion
}