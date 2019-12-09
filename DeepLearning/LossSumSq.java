package DeepLearning;

//region Imports


//endregion

public class LossSumSq implements ILossFunction
{

	//region ILossFunction Implementation
	public final double Measure(Neuron exp, Neuron act)
	{
		double sum = 0;
		for (int i = 0; i < exp.W.getLength(); i++)
		{
			double errDelta = act.W[i] - exp.W[i];
			sum += 0.5 * errDelta * errDelta;
		}
		return sum;
	}

	public final void Backward(Neuron exp, Neuron act)
	{
		for (int i = 0; i < exp.W.getLength(); i++)
		{
			double errDelta = act.W[i] - exp.W[i];
			act.Dw[i] += errDelta;
		}
	}

	//endregion
}