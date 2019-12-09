package DeepLearning;

import java.util.*;

//region Imports


//endregion


public final class NNFactory
{

	//region Public Methods

	public static NN CreateFForward(int inpDim, int hidDim, int outDim, int hidLayers, ActFunction hidActFun, ActFunction outActFun, double stDev)
	{
		ILayer hid0Layer = new LayerFF(inpDim, hidDim, stDev, hidActFun);
		ILayer hid1Layer = new LayerFF(hidDim, hidDim, stDev, hidActFun);
		ILayer outLayer = new LayerFF(hidDim, outDim, stDev, outActFun);
		ArrayList<ILayer> layers = new ArrayList<ILayer>();
		return CreateNN(hidLayers, null, hid0Layer, hid1Layer, outLayer);
	}

	public static NN CreateRnn(int inpDim, int hidDim, int outDim, int hidLayers, ActFunction hidActFun, ActFunction octActFun, double stDev)
	{
		ILayer hid0Layer = new LayerRNN(inpDim, hidDim, stDev, hidActFun);
		ILayer hid1Layer = new LayerRNN(hidDim, hidDim, stDev, hidActFun);
		ILayer outLayer = new LayerFF(hidDim, outDim, stDev, octActFun);
		ArrayList<ILayer> layers = new ArrayList<ILayer>();
		return CreateNN(hidLayers, null, hid0Layer, hid1Layer, outLayer);
	}

	public static NN CreateGru(int inpDim, int hidDim, int outDim, int hidLayers, ActFunction outActFun, double stDev)
	{
		ILayer hid0Layer = new LayerGRU(inpDim, hidDim, stDev);
		ILayer hid1Layer = new LayerGRU(hidDim, hidDim, stDev);
		ILayer outLayer = new LayerFF(hidDim, outDim, stDev, outActFun);
		ArrayList<ILayer> layers = new ArrayList<ILayer>();
		return CreateNN(hidLayers, null, hid0Layer, hid1Layer, outLayer);
	}

	public static NN CreateLstm(int inpDim, int hidDim, int outDim, int hidLayers, ActFunction outActFun, double stDev)
	{
		ILayer hid0Layer = new LayerLSTM(inpDim, hidDim, stDev);
		ILayer hid1Layer = new LayerLSTM(hidDim, hidDim, stDev);
		ILayer outLayer = new LayerFF(hidDim, outDim, stDev, outActFun);
		ArrayList < ILayer> layers = new ArrayList<ILayer>();
		return CreateNN(hidLayers, null, hid0Layer, hid1Layer, outLayer);
	}

	public static NN CreateLstmBneck(int inpDim, int hidDim, int outDim, int hidLayers, int bottleneckDim, ActFunction outActFun, double stDev)
	{ //with input bottleneck
		ILayer inpLayer = new LayerLin(inpDim, bottleneckDim, stDev);
		ILayer hid0Layer = new LayerLSTM(bottleneckDim, hidDim, stDev);
		ILayer hid1Layer = new LayerLSTM(hidDim, hidDim, stDev);
		ILayer outLayer = new LayerFF(hidDim, outDim, stDev, outActFun);
		return CreateNN(hidLayers, inpLayer, hid0Layer, hid1Layer, outLayer);
	}

	//endregion

	//region Private Methods

	private static NN CreateNN(int hidLayers, ILayer inpLayer, ILayer hid0Layer, ILayer hid1Layer, ILayer outLayer)
	{
		ArrayList<ILayer> layers = new ArrayList<ILayer>();
		if (inpLayer != null)
		{
			layers.add(inpLayer);
		}
		for (int h = 0; h < hidLayers; h++)
		{
			if (h == 0)
			{
				layers.add(hid0Layer);
			}
			else
			{
				layers.add(hid1Layer);
			}
		}
		layers.add(outLayer);
		return new NN(layers);
	}

	//endregion
}