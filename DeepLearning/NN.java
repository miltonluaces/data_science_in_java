package DeepLearning;

import java.util.*;

//endregion

public class NN implements INN, Serializable
{

	//region Fields

	private static long uid = 1;
	private ArrayList<ILayer> layers;

	//endregion

	//region Constructor
	public NN(ArrayList<ILayer> layers)
	{
		this.layers = layers;
	}

	//endregion

	//region INetwork Implementation

	public final Neuron Activate(Neuron input, NNActivation nnAct)
	{
		Neuron prev = input;
		for (ILayer layer : layers)
		{
			prev = layer.Activate(prev, nnAct);
		}
		return prev;
	}

	public final ArrayList<Neuron> GetParams()
	{
		ArrayList<Neuron> pars = new ArrayList<Neuron>();
		for (ILayer layer : layers)
		{
			pars.addAll(layer.GetParams());
		}
		return pars;
	}

	public final void Reset()
	{
		for (ILayer layer : layers)
		{
			layer.Reset();
		}
	}

	//endregion

	//region Persistence

	public static NN Read(String filePath)
	{
		if ((new java.io.File(filePath)).isFile())
		{
			try (FileStream fs = new FileStream(filePath, FileMode.Open))
			{
				BinaryFormatter formatter = new BinaryFormatter();
				formatter.AssemblyFormat = System.Runtime.Serialization.Formatters.FormatterAssemblyStyle.Simple;
				if (fs.getLength() == 0)
				{
					return null;
				}
				return ((NN)formatter.Deserialize(fs));
			}
		}
		throw new FileNotFoundException(filePath);
	}

	public static void Write(NN nn, String filePath)
	{
		try (Stream stream = File.Open(filePath, FileMode.Create))
		{
			BinaryFormatter bformatter = new BinaryFormatter();
			bformatter.AssemblyFormat = System.Runtime.Serialization.Formatters.FormatterAssemblyStyle.Simple;
			bformatter.Serialize(stream, nn);
		}
	}

	//endregion
}