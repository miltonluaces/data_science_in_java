package DeepLearning;

import java.util.*;

public class Data
{

	private Neuron Input;
	public final Neuron getInput()
	{
		return Input;
	}
	public final void setInput(Neuron value)
	{
		Input = value;
	}
	private Neuron ExpOutput;
	public final Neuron getExpOutput()
	{
		return ExpOutput;
	}
	public final void setExpOutput(Neuron value)
	{
		ExpOutput = value;
	}

	public Data()
	{
	}
	public Data(double[] inp, double[] expOut)
	{
		this.setInput(new Neuron(inp));
		if (expOut != null)
		{
			this.setExpOutput(new Neuron(expOut));
		}
	}

	@Override
	public String toString()
	{
		String str = "";
		for (int i = 0; i < getInput().W.length; i++)
		{
//C# TO JAVA CONVERTER TODO TASK: The '0:N5' format specifier is not converted to Java:
			str += String.format("{0:N5}", getInput().W[i]) + "\t";
		};
		str += "\t->\t";
		if (getExpOutput() != null)
		{
			for (int i = 0; i < getExpOutput().W.length; i++)
			{
//C# TO JAVA CONVERTER TODO TASK: The '0:N5' format specifier is not converted to Java:
				str += String.format("{0:N5}", getExpOutput().W[i]) + "\t";
			}
		}
		else
		{
			str += "   \t";
		}
		return str;
	}
}