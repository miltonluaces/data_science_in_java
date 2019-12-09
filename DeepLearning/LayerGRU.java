package DeepLearning;

import java.util.*;

//region Imports


//endregion


 public class LayerGRU implements ILayer, Serializable
 {

	//region Fields

	private static long uid = 1L;

	private int inpDim;
	private int outDim;

	private Neuron hmix;
	private Neuron hHmix;
	private Neuron bmix;
	private Neuron hnew;
	private Neuron hHnew;
	private Neuron bnew;
	private Neuron hreset;
	private Neuron hHreset;
	private Neuron breset;

	private Neuron context;

	private ActFunction fMix;
	private ActFunction fReset;
	private ActFunction fNew;

	//endregion

	//region Constructor

	public LayerGRU(int inputDim, int outputDim, double stDev)
	{

		this.inpDim = inputDim;
		this.outDim = outputDim;

		fMix = new ActSigm();
		fReset = new ActSigm();
		fNew = new ActTanHyp();

		hmix = Neuron.Random(outputDim, inputDim, stDev);
		hHmix = Neuron.Random(outputDim, outputDim, stDev);
		bmix = new Neuron(outputDim);
		hnew = Neuron.Random(outputDim, inputDim, stDev);
		hHnew = Neuron.Random(outputDim, outputDim, stDev);
		bnew = new Neuron(outputDim);
		hreset = Neuron.Random(outputDim, inputDim, stDev);
		hHreset = Neuron.Random(outputDim, outputDim, stDev);
		breset = new Neuron(outputDim);
	}

	//endregion

	//region ILayer Implementation

	public final Neuron Activate(Neuron input, NNActivation nnAct)
	{

	Neuron sum0 = nnAct.Prod(hmix, input);
	Neuron sum1 = nnAct.Prod(hHmix, context);
	Neuron actMix = nnAct.Activate(fMix, nnAct.Add(nnAct.Add(sum0, sum1), bmix));

	Neuron sum2 = nnAct.Prod(hreset, input);
	Neuron sum3 = nnAct.Prod(hHreset, context);
	Neuron actReset = nnAct.Activate(fReset, nnAct.Add(nnAct.Add(sum2, sum3), breset));

	Neuron sum4 = nnAct.Prod(hnew, input);
	Neuron gatedContext = nnAct.Elmul(actReset, context);
	Neuron sum5 = nnAct.Prod(hHnew, gatedContext);
	Neuron actNewPlusGatedContext = nnAct.Activate(fNew, nnAct.Add(nnAct.Add(sum4, sum5), bnew));

	Neuron memvals = nnAct.Elmul(actMix, context);
	Neuron newvals = nnAct.Elmul(nnAct.OneMinus(actMix), actNewPlusGatedContext);
	Neuron output = nnAct.Add(memvals, newvals);

	context = output; //rollover activations for next iteration
	return output;
	}

public final ArrayList<Neuron> GetParams()
{
	Neuron[] parsArray = {hmix, hHmix, bmix, hnew, hHnew, bnew, hreset, hHreset, breset};
	return new ArrayList<Neuron>(parsArray);
}

public final void Reset()
{
	context = new Neuron(outDim);
}

//endregion

 }