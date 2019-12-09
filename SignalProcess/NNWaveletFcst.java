package SignalProcess;

import NeuralNetworks.*;
import NumCalc.*;
import ClassicalStat.*;
import java.util.*;

//region Imports


//endregion


/**  Class for Wavelet neural network (mixed technique) forecasting 
*/
public class NNWaveletFcst
{

	//region Fields

	private Wavelets ws;
	private TDNN tdnn;
	private ArrayList<ArrayList<Double>> fcsts;
	private ArrayList<Double> fcstTot;
	private ArrayList<Double> leadTimeFcsts;
	private ArrayList<Double> fcstability;
	private int fcstHorizon;

	//endregion

	//region Constructor

	/**  Constructor 
	*/
	public NNWaveletFcst()
	{
		this.ws = new Wavelets();
		this.fcsts = new ArrayList<ArrayList<Double>>();
		this.leadTimeFcsts = new ArrayList<Double>();
		this.fcstability = new ArrayList<Double>();
		this.fcstTot = new ArrayList<Double>();
	}

	//endregion

	//region Properties

	/**  Original time series 
	*/
	public final ArrayList<Double> getOrigData()
	{
		return ws.OrigData;
	}

	/**  Current time series 
	*/
	public final ArrayList<Double> getData()
	{
		return ws.Data;
	}

	/**  List of wavelets of the model 
	*/
	public final Wavelets getWs()
	{
		return ws;
	}

	/**  Time delayed neural network for wnn calculation 
	*/
	public final TDNN getTdnn()
	{
		return tdnn;
	}

	/**  Lead time forecast time series 
	*/
	public final ArrayList<Double> getLeadTimeFcsts()
	{
		return leadTimeFcsts;
	}

	/**  Forecasting horizon 
	*/
	public final int getFcstHorizon()
	{
		return fcstHorizon;
	}

	//endregion     

	//region Forecast Getters

	/**  Get forecast of a particular spectrum 
	 @param spec spectrum index 
	 @return  Forecast time series 
	*/
	public final ArrayList<Double> GetFcst(int spec)
	{
		if (spec >= fcsts.size())
		{
			throw new RuntimeException(Strings.Spec_must_be_less_or_equal_than_count);
		}
		return fcsts.get(spec);
	}

	/**  Get forecast of all spectrums until one 
	 @param spec last spectrum index 
	 @return  Forecast time series 
	*/
	public final ArrayList<Double> GetTotFcstUntil(int spec)
	{
		if (spec > fcsts.size())
		{
			throw new RuntimeException(Strings.Spec_must_be_less_or_equal_than_count);
		}
		ArrayList<Double> untilFcst = new ArrayList<Double>();
		for (int i = 0;i < fcstHorizon;i++)
		{
			untilFcst.add(0.00);
		}
		for (int i = 0;i < spec;i++)
		{
			for (int j = 0;j < fcstHorizon;j++)
			{
				untilFcst.set(j, untilFcst.get(j) + fcsts.get(i).get(j));
			}
		}
		return untilFcst;
	}

	/**  Get total forecast 
	 @return  forecast time series 
	*/
	public final ArrayList<Double> GetTotFcst()
	{
		return GetTotFcstUntil(fcsts.size());
	}

	/**  Get forecastability level for a particular spectrum 
	 @param spec spectrum index 
	 @return  forecastability level 
	*/
	public final double GetFcstability(int spec)
	{
		return fcstability.get(spec - 1);
	}

	//endregion

	//region Public Methods

	/**  Load data for calculation 
	 @param data data time series 
	 @param group number of elements per group (previous grouping) 
	*/
	public final void LoadData(ArrayList<Double> data, int group, boolean continuous)
	{
		ws.Group(data, 0, group, continuous);
		ws.CalcFreqDescomposition();
	}

	/**  Main calculate method 
	 @param fcstHorizon forecasting horizon 
	 @param epochs number of epochs 
	 @param serviceLevel percentage of service level 
	 @param errorWindow width of the error window 
	*/
	public final void Calculate(int fcstHorizon, int epochs, double serviceLevel, int errorWindow)
	{
		this.fcstHorizon = fcstHorizon;
		fcsts.clear();
		ArrayList<Double> fcst = new ArrayList<Double>();
		int mw = 1;
		int maxMw = 16;
		ArrayList<Double> desc = null;
		double min, max;
		fcst = new ArrayList<Double>();
		desc = ws.GetFreqDiffDescomposition(0);
		for (int i = 0;i < fcstHorizon;i++)
		{
			fcst.add(desc.get(desc.size() - 1));
		}
		fcsts.add(fcst);

		//forecast
		for (int s = 1;s < ws.Coeffs.size();s++)
		{
			fcst = new ArrayList<Double>();
			if (mw < maxMw)
			{
				mw += 2;
			}
			desc = ws.GetFreqDiffDescomposition(s);
			min = Double.MAX_VALUE;
			max = -Double.MAX_VALUE;
			for (double val : desc)
			{
				if (val < min)
				{
					min = val;
				}
				if (val > max)
				{
					max = val;
				}
			}
			if (min == max)
			{
				for (int i = 0;i < fcstHorizon;i++)
				{
					fcst.add(min);
				}
				fcsts.add(fcst);
				continue;
			}
			if (min < 0)
			{
				for (int i = 0;i < desc.size();i++)
				{
					desc.set(i, desc.get(i) - min);
				}
			}
			tdnn = new TDNN(mw);
			tdnn.LoadSerie(desc, 0);
			double[] fcstArr = tdnn.Process(epochs, fcstHorizon, errorWindow);
			fcst = new ArrayList<Double>();
			if (min >= 0)
			{
				for (double val : fcstArr)
				{
					fcst.add(val);
				}
			}
			else
			{
				for (double val : fcstArr)
				{
					fcst.add(val + min);
				}
			}
			fcsts.add(fcst);
		}
		fcstTot.clear();
		for (int i = 0;i < fcstHorizon;i++)
		{
			fcstTot.add(0.0);
			for (ArrayList<Double> fcsT : fcsts)
			{
				fcstTot.set(i, fcstTot.get(i) + fcsT.get(i));
			}
		}

		//errors 
		double plus, frcst;
		if (tdnn == null || tdnn.Errors == null)
		{
			plus = 0.0;
			leadTimeFcsts.clear();
			for (int i = 0;i < fcstHorizon;i++)
			{
				frcst = fcstTot.get(i) * (1.0 + plus);
				if (frcst < 0)
				{
					frcst = 0.0;
				}
				leadTimeFcsts.add(frcst);
			}
			return;
		}
		DensityEstim dEst = new DensityEstim(1, 30, 100);
		ArrayList<Double> errors = tdnn.Errors;
		Histogram hist = new Histogram(100);
		hist.LoadData(errors);
		dEst.LoadDist(hist.GetFreqs());
		dEst.SetMaxInt();
		double slFixed = serviceLevel + 4.00;
		if (slFixed > 99.0)
		{
			slFixed = 99.0;
		}
		double error = dEst.CalculatePercentile(slFixed);


		plus = 0.0;
		leadTimeFcsts.clear();
		for (int i = 0;i < fcstHorizon;i++)
		{
			frcst = (fcstTot.get(i) + error * (1.0 + plus));
			if (frcst < 0)
			{
				frcst = 0.0;
			}
			leadTimeFcsts.add(frcst);
		}
	}

	//endregion
}