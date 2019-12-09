package DeepLearning;

//region Imports


//endregion

public class Entropy implements ILossFunction
{

	//region ILossFunction Implementation

	public final double Measure(Neuron exp, Neuron act)
	{
		double crossEntropy = 0.0;
		for (int i = 0; i < act.W.getLength(); i++)
		{
			crossEntropy -= (exp.W[i] * Math.log(act.W[i] + 1e-15)) + ((1 - exp.W[i]) * Math.log((1 + 1e-15) - act.W[i]));
		}
		return crossEntropy;
	}

	public final void Backward(Neuron exp, Neuron act)
	{
		throw new RuntimeException("not implemented");
	}

	//endregion
}