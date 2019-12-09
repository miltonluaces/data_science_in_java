package DeepLearning;

import java.util.*;

//region Imports


//endregion

public class LossSoftmax implements ILossFunction
{

	//region ILossFunction Implementation

	public final double Measure(Neuron exp, Neuron logprobs)
	{
		int targetIndex = GetExpOutIndex(exp);
		Neuron probs = GetSoftmaxProbs(logprobs, 1.0);
		double cost = -Math.log(probs.W[targetIndex]);
		return cost;
	}

	public final void Backward(Neuron exp, Neuron logprobs)
	{
		int targetIndex = GetExpOutIndex(exp);
		Neuron probs = GetSoftmaxProbs(logprobs, 1.0);
		for (int i = 0; i < probs.W.getLength(); i++)
		{
			logprobs.Dw[i] = probs.W[i];
		}
		logprobs.Dw[targetIndex] -= 1;
	}

	//endregion

	//region Public Static Methods

	public static double CalculateMedianPerplexity(ILayer layer, ArrayList<ArrayList<Data>> seqs)
	{
		double temperature = 1.0;
		ArrayList<Double> ppls = new ArrayList<Double>();
		for (ArrayList<Data> seq : seqs)
		{
			double n = 0;
			double neglog2Ppl = 0;

			NNActivation nnAct = new NNActivation(false);
			layer.Reset();
			for (Data step : seq)
			{
				Neuron logprobs = layer.Activate(step.Input, nnAct);
				Neuron probs = GetSoftmaxProbs(logprobs, temperature);
				int targetIndex = GetExpOutIndex(step.ExpOutput);
				double probOfCorrect = probs.W[targetIndex];
				double log2Prob = Math.log(probOfCorrect) / Math.log(2); //change-of-base
				neglog2Ppl += -log2Prob;
				n += 1;
			}

			n -= 1; //don't count first symbol of sentence
			double ppl = Math.pow(2, (neglog2Ppl / (n - 1)));
			ppls.add(ppl);
		}
		return Median(ppls);
	}

	public static Neuron GetSoftmaxProbs(Neuron logprobs, double temperature)
	{
		Neuron probs = new Neuron(logprobs.W.getLength());
		if (temperature != 1.0)
		{
			for (int i = 0; i < logprobs.W.getLength(); i++)
			{
				logprobs.W[i] /= temperature;
			}
		}
		double maxval = Double.NEGATIVE_INFINITY;
		for (int i = 0; i < logprobs.W.getLength(); i++)
		{
			if (logprobs.W[i] > maxval)
			{
				maxval = logprobs.W[i];
			}
		}
		double sum = 0;
		for (int i = 0; i < logprobs.W.getLength(); i++)
		{
			probs.W[i] = Math.exp(logprobs.W[i] - maxval); //all inputs to exp() are non-positive
			sum += probs.W[i];
		}
		for (int i = 0; i < probs.W.getLength(); i++)
		{
			probs.W[i] /= sum;
		}
		return probs;
	}

	//endregion

	//region Private Auxiliar Methods
	public static double Median(ArrayList<Double> vals)
	{
//C# TO JAVA CONVERTER TODO TASK: There is no Java equivalent to LINQ queries:
		vals = vals.OrderBy(a -> a).ToList();
		int mid = vals.size() / 2;
		if (vals.size() % 2 == 1)
		{
			return vals.get(mid);
		}
		else
		{
			return (vals.get(mid - 1) + vals.get(mid)) / 2;
		}
	}

	private static int GetExpOutIndex(Neuron expOut)
	{
		for (int i = 0; i < expOut.W.getLength(); i++)
		{
			if (expOut.W[i] == 1.0)
			{
				return i;
			}
		}
		throw new RuntimeException("no expected output index selected");
	}

	//endregion
}