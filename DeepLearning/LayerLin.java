package DeepLearning;

import java.util.*;

//region Imports


//endregion


public class LayerLin implements ILayer, Serializable
{

	//region Fields

	private static long uid = 1L;

	private Neuron W; //no biases

	//endregion

	//region Constructor

	public LayerLin(int inpDim, int outDim, double stDev)
	{
		W = Neuron.Random(outDim, inpDim, stDev);
	}

	//endregion

	//region ILayer Implementation

	public final Neuron Activate(Neuron input, NNActivation nnAct)
	{
		return nnAct.Prod(W, input);
	}

	public final ArrayList<Neuron> GetParams()
	{
		ArrayList<Neuron> pars = new ArrayList<Neuron>();
		pars.add(W);
		return pars;
	}

	public final void Reset()
	{
		//Not implemented
	}

	//endregion
}