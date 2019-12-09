package DeepLearning;

import java.util.*;

//region Imports


//endregion


public class LayerLSTM implements ILayer, Serializable
{

	//region Fields

	private static long uid = 1L;

	private int inpDim;
	private int outDim;

	private Neuron inpWx;
	private Neuron inpWh;
	private Neuron inpBias;
	private Neuron forgWx;
	private Neuron forgWh;
	private Neuron forgBias;
	private Neuron outWx;
	private Neuron outWh;
	private Neuron outBias;
	private Neuron cellWx;
	private Neuron cellWh;
	private Neuron cellBias;

	private Neuron hiddenContext;
	private Neuron cellContext;

	private ActFunction inputActFunction = new ActSigm();
	private ActFunction forgetActFunction = new ActSigm();
	private ActFunction outputActFunction = new ActSigm();
	private ActFunction cellInpActFunction = new ActTanHyp();
	private ActFunction cellOutActFunction = new ActTanHyp();

	//endregion

	//region Constructor
	public LayerLSTM(int inpDim, int outDim, double stDev)
	{
		this.inpDim = inpDim;
		this.outDim = outDim;

		this.inpWx = Neuron.Random(outDim, inpDim, stDev);
		this.inpWh = Neuron.Random(outDim, outDim, stDev);
		this.inpBias = new Neuron(outDim);

		this.forgWx = Neuron.Random(outDim, inpDim, stDev);
		this.forgWh = Neuron.Random(outDim, outDim, stDev);
		this.forgBias = Neuron.Rep(outDim, 1, 1.0); //set forget bias to 1.0, as described here: http://jmlr.org/proceedings/papers/v37/jozefowicz15.pdf

		this.outWx = Neuron.Random(outDim, inpDim, stDev);
		this.outWh = Neuron.Random(outDim, outDim, stDev);
		this.outBias = new Neuron(outDim);

		this.cellWx = Neuron.Random(outDim, inpDim, stDev);
		this.cellWh = Neuron.Random(outDim, outDim, stDev);
		this.cellBias = new Neuron(outDim);
	}

	//endregion

	//region ILayer Implementation

	public final Neuron Activate(Neuron input, NNActivation nnAct)
	{

		//input gate
		Neuron sum0 = nnAct.Prod(inpWx, input);
		Neuron sum1 = nnAct.Prod(inpWh, hiddenContext);
		Neuron inputGate = nnAct.Activate(inputActFunction, nnAct.Add(nnAct.Add(sum0, sum1), inpBias));

		//forget gate
		Neuron sum2 = nnAct.Prod(forgWx, input);
		Neuron sum3 = nnAct.Prod(forgWh, hiddenContext);
		Neuron forgetGate = nnAct.Activate(forgetActFunction, nnAct.Add(nnAct.Add(sum2, sum3), forgBias));

		//output gate
		Neuron sum4 = nnAct.Prod(outWx, input);
		Neuron sum5 = nnAct.Prod(outWh, hiddenContext);
		Neuron outputGate = nnAct.Activate(outputActFunction, nnAct.Add(nnAct.Add(sum4, sum5), outBias));

		//write operation on cells
		Neuron sum6 = nnAct.Prod(cellWx, input);
		Neuron sum7 = nnAct.Prod(cellWh, hiddenContext);
		Neuron cellInput = nnAct.Activate(cellInpActFunction, nnAct.Add(nnAct.Add(sum6, sum7), cellBias));

		//compute new cell activation
		Neuron retainCell = nnAct.Elmul(forgetGate, cellContext);
		Neuron writeCell = nnAct.Elmul(inputGate, cellInput);
		Neuron cellAct = nnAct.Add(retainCell, writeCell);

		//compute hidden state as gated, saturated cell activations
		Neuron output = nnAct.Elmul(outputGate, nnAct.Activate(cellOutActFunction, cellAct));

		//rollover activations for next iteration
		hiddenContext = output;
		cellContext = cellAct;

		return output;
	}

	public final ArrayList<Neuron> GetParams()
	{
		Neuron[] parsArray = {inpWx, inpWh, inpBias, forgWx, forgWh, forgBias, outWx, outWh, outBias, cellWx, cellWh, cellBias};
		return new ArrayList<Neuron>(parsArray);
	}

	public final void Reset()
	{
		hiddenContext = new Neuron(outDim);
		cellContext = new Neuron(outDim);
	}

	//endregion
}