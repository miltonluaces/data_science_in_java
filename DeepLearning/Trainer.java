package DeepLearning;

import java.util.*;

//region Imports


//endregion

public class Trainer
{

	//region Fields

	public static double DecayRate = 0.999;
	public static double SmoothEpsilon = 1e-8;
	public static double GradientClipValue = 5;
	public static double Regularization = 0.000001; // L2 regularization strength

	//endregion

	//region Public Methods

	public static double Train(int epochs, double learningRate, INN network, Dataset data, int reportFreq)
	{
		return Train(epochs, learningRate, network, data, reportFreq, false, false, null);
	}

	public static double Train(int epochs, double learningRate, INN network, Dataset ds, int reportFreq, boolean initFromSaved, boolean overwriteSaved, String savePath)
	{
		System.out.println("--------------------------------------------------------------");
		if (initFromSaved)
		{
			System.out.println("initializing network from saved state...");
			try
			{
				network = NN.Read(savePath);
			}
			catch (RuntimeException e)
			{
				System.out.println("Error. Unable to load from a saved state. WARNING: " + e.getMessage());
			}
		}
		double result = 1.0;
		for (int epoch = 0; epoch < epochs; epoch++)
		{
			String show = "epoch[" + (epoch + 1) + "/" + epochs + "]";

			double reportedLossTrain = Pass(learningRate, network, ds.getTrainSet(), true, ds.getLossTrain(), ds.getLossReport());
			result = reportedLossTrain;
			if (Double.isNaN(reportedLossTrain) || Double.isInfinite(reportedLossTrain))
			{
				throw new RuntimeException("WARNING: invalid value for training loss. Try lowering learning rate.");
			}
			double reportedLossValidation = 0;
			double reportedLossTesting = 0;
			if (ds.getValidSet() != null)
			{
				reportedLossValidation = Pass(learningRate, network, ds.getValidSet(), false, ds.getLossTrain(), ds.getLossReport());
				result = reportedLossValidation;
			}
			if (ds.getTestSet() != null)
			{
				reportedLossTesting = Pass(learningRate, network, ds.getTestSet(), false, ds.getLossTrain(), ds.getLossReport());
				result = reportedLossTesting;
			}
//C# TO JAVA CONVERTER TODO TASK: The '0:N5' format specifier is not converted to Java:
			show += "\ttrain loss = " + String.format("{0:N5}", reportedLossTrain);
			if (ds.getValidSet() != null)
			{
//C# TO JAVA CONVERTER TODO TASK: The '0:N5' format specifier is not converted to Java:
				show += "\tvalid loss = " + String.format("{0:N5}", reportedLossValidation);
			}
			if (ds.getTestSet() != null)
			{
//C# TO JAVA CONVERTER TODO TASK: The '0:N5' format specifier is not converted to Java:
				show += "\ttest loss  = " + String.format("{0:N5}", reportedLossTesting);
			}
			System.out.println(show);

			if (overwriteSaved)
			{
				NN.Write((NN)network, savePath);
			}

			if (reportedLossTrain == 0 && reportedLossValidation == 0)
			{
				System.out.println("--------------------------------------------------------------");
				System.out.println("\nDONE.");
				break;
			}
		}
		return result;
	}

	//endregion

	//region Private Methods
	private static double Pass(double learningRate, INN network, ArrayList<ArrayList<Data>> sequences, boolean applyTraining, ILossFunction lossTraining, ILossFunction lossReporting)
	{
		double num = 0;
		double den = 0;

		for (ArrayList<Data> sequence : sequences)
		{
			network.Reset();
			NNActivation graph = new NNActivation(applyTraining);
			for (Data step : sequence)
			{
				Neuron output = network.Activate(step.getInput(), graph);
				if (step.getExpOutput() != null)
				{
					double loss = lossReporting.Measure(output, step.getExpOutput());
					if (Double.isNaN(loss) || Double.isInfinite(loss))
					{
						return loss;
					}
					num += loss;
					den++;
					if (applyTraining)
					{
						lossTraining.Backward(output, step.getExpOutput());
					}
				}
			}
			ArrayList<ArrayList<Data>> thisSequence = new ArrayList<ArrayList<Data>>();
			thisSequence.add(sequence);
			if (applyTraining)
			{
				graph.Backward(); //backprop dw values
				UpdateParams(network, learningRate); //update params
			}
		}
		return num / den;
	}

	private static void UpdateParams(INN network, double stepSize)
	{
		for (Neuron m : network.GetParams())
		{
			for (int i = 0; i < m.W.length; i++)
			{

				// rmsprop adaptive learning rate
				double mdwi = m.Dw[i];
				m.Cache[i] = m.Cache[i] * DecayRate + (1 - DecayRate) * mdwi * mdwi;

				// gradient clip
				if (mdwi > GradientClipValue)
				{
					mdwi = GradientClipValue;
				}
				if (mdwi < -GradientClipValue)
				{
					mdwi = -GradientClipValue;
				}

				// update (and regularize)
				m.W[i] += -stepSize * mdwi / Math.sqrt(m.Cache[i] + SmoothEpsilon) - Regularization * m.W[i];
				m.Dw[i] = 0;
			}
		}
	}

	//endregion

}