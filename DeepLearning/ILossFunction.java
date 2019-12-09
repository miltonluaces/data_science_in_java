package DeepLearning;

//region Imports


//endregion


public interface ILossFunction
{
	double Measure(Neuron exp, Neuron act);
	void Backward(Neuron exp, Neuron act);
}