package DeepLearning;

import java.util.*;

//region Imports


//endregion


//region Interface INetwork

public interface INN
{

	Neuron Activate(Neuron input, NNActivation graph);
	ArrayList<Neuron> GetParams();
	void Reset();

}