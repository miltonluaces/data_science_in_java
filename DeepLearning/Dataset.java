package DeepLearning;

import java.util.*;

//region Imports


//endregion


public class Dataset
{

	public Dataset()
	{
	}
	private int InpDim;
	public final int getInpDim()
	{
		return InpDim;
	}
	public final void setInpDim(int value)
	{
		InpDim = value;
	}
	private int OutDim;
	public final int getOutDim()
	{
		return OutDim;
	}
	public final void setOutDim(int value)
	{
		OutDim = value;
	}
	private ActFunction ActFunction;
	public final ActFunction getActFunction()
	{
		return ActFunction;
	}
	public final void setActFunction(ActFunction value)
	{
		ActFunction = value;
	}
	private ILossFunction LossTrain;
	public final ILossFunction getLossTrain()
	{
		return LossTrain;
	}
	public final void setLossTrain(ILossFunction value)
	{
		LossTrain = value;
	}
	private ILossFunction LossReport;
	public final ILossFunction getLossReport()
	{
		return LossReport;
	}
	public final void setLossReport(ILossFunction value)
	{
		LossReport = value;
	}
	private ArrayList<ArrayList<Data>> TrainSet;
	public final ArrayList<ArrayList<Data>> getTrainSet()
	{
		return TrainSet;
	}
	public final void setTrainSet(ArrayList<ArrayList<Data>> value)
	{
		TrainSet = value;
	}
	private ArrayList<ArrayList<Data>> ValidSet;
	public final ArrayList<ArrayList<Data>> getValidSet()
	{
		return ValidSet;
	}
	public final void setValidSet(ArrayList<ArrayList<Data>> value)
	{
		ValidSet = value;
	}
	private ArrayList<ArrayList<Data>> TestSet;
	public final ArrayList<ArrayList<Data>> getTestSet()
	{
		return TestSet;
	}
	public final void setTestSet(ArrayList<ArrayList<Data>> value)
	{
		TestSet = value;
	}
}