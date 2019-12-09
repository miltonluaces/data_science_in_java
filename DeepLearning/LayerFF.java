package DeepLearning;

import java.util.*;

//region Imports


//endregion


public class LayerFF implements ILayer, Serializable
{

	//region Fields

	private static long uid = 1;

	private Neuron w;
	private Neuron b;
	private ActFunction neu;

	//endregion

	//region Constructor

	public LayerFF(int inpDim, int outDim, double stDev, ActFunction actFun)
	{
		w = Neuron.Random(outDim, inpDim, stDev);
		b = new Neuron(outDim);
		this.neu = actFun;
	}

	//endregion

	//region ILayer Implementation

	public final Neuron Activate(Neuron input, NNActivation nnAct)
	{
		Neuron sum = nnAct.Add(nnAct.Prod(w, input), b);
		Neuron res = nnAct.Activate(neu, sum);
		return res;
	}

	public final ArrayList<Neuron> GetParams()
	{
		ArrayList<Neuron> pars = new ArrayList<Neuron>();
		pars.add(w);
		pars.add(b);
		return pars;
	}

	public final void Reset()
	{
		//Not implemented
	}

	//endregion
}