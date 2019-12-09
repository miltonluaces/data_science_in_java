package DeepLearning;

import java.util.*;

//region Imports


//endregion


public interface ILayer
{

	Neuron Activate(Neuron input, NNActivation nnAct);
	ArrayList<Neuron> GetParams();
	void Reset();
}