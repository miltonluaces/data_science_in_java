package DeepLearning;

//region Imports


//endregion

public class LossArgMax implements ILossFunction
{

	//region ILossFunction Implementation

	public final double Measure(Neuron exp, Neuron act)
	{
		if (act.W.getLength() != exp.W.getLength())
		{
			throw new RuntimeException("mismatch");
		}

		double maxAct = Double.POSITIVE_INFINITY;
		double maxExp = Double.NEGATIVE_INFINITY;
		int iMaxAct = -1;
		int iMaxExp = -1;
		for (int i = 0; i < act.W.getLength(); i++)
		{
			if (act.W[i] > maxAct)
			{
				maxAct = act.W[i];
				iMaxAct = i;
			}
			if (exp.W[i] > maxExp)
			{
				maxExp = exp.W[i];
				iMaxExp = i;
			}
		}
		return (iMaxAct == iMaxExp)? 0 : 1;
	}

	public final void Backward(Neuron exp, Neuron act)
	{
		throw new RuntimeException("not implemented");
	}

	//endregion
}