package DeepLearning;

import java.util.*;

//region Imports


//endregion


public class LayerRNN implements ILayer, Serializable
{

	//region Fields

	private static long uid = 1;

	private int inputDim;
	private int outputDim;

	private Neuron W;
	private Neuron B;

	private Neuron context;

	private ActFunction actFun;

	//endregion

	//region Constructor
	public LayerRNN(int inpDim, int outDim, double stDev, ActFunction hidActFun)
	{
		this.inputDim = inpDim;
		this.outputDim = outDim;
		this.actFun = hidActFun;
		this.W = Neuron.Random(outDim, inpDim + outDim, stDev);
		this.B = new Neuron(outDim);
	}

	//endregion

	//region ILayer Implementation

	public final Neuron Activate(Neuron input, NNActivation nnAct)
	{
		Neuron concat = nnAct.Concat(input, context);
		Neuron sum = nnAct.Prod(W, concat);
		sum = nnAct.Add(sum, B);
		Neuron output = nnAct.Activate(actFun, sum); //rollover activations for next iteration
		context = output;
		return output;
	}

	public final ArrayList<Neuron> GetParams()
	{
		Neuron[] parsArray = {W, B};
		return new ArrayList<Neuron>(parsArray);
	}

	public final void Reset()
	{
		context = new Neuron(outputDim);
	}

	//endregion
}