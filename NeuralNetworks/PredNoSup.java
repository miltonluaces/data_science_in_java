package NeuralNetworks;

import java.util.*;

//region Imports


//endregion

/**  Class non-supervised prediction 
*/
public class PredNoSup
{

	//region Fields

	private TDNN tdnn;
	private int mw;
	private int ciclos;
	private int fcstHorizon;
	private double umbral;
	private int nUltPeriodos;
	private double fc;

	//endregion

	//region Constructor

	/**  Constructor 
	 @param mw moving window 
	 @param ciclos epochs 
	 @param fcstHorizon forecasting horizon 
	 @param umbral threshold used in calculation 
	 @param nUltPeriodos number of last periods 
	 @param fc load factor (proportion of zeros) 
	*/
	public PredNoSup(int mw, int ciclos, int fcstHorizon, double umbral, int nUltPeriodos, double fc)
	{
		this.mw = mw;
		this.ciclos = ciclos;
		this.fcstHorizon = fcstHorizon;
		this.umbral = umbral;
		this.nUltPeriodos = nUltPeriodos;
		this.fc = fc;
	}

	//endregion

	//region Public Methods

	/**  Devuelve si la serie es o no forecastable 
	 @param serie the time series 
	 @return  if it is forecastable 
	*/
	public final boolean IsForecastable(ArrayList<Double> serie)
	{
		if (serie.size() < mw * 2)
		{
			return false;
		}
		if (GetFc(serie) < fc)
		{
			return false;
		}
		if (HasUltPerZero(serie))
		{
			return false;
		}
		/*Console.WriteLine("Determinado por RNA"); */ return HowForecastable(serie) >= umbral;
	}

	/**  Devuelve un valor entre 0 y 1 que representa el nivel de pronosticabilidad 
	 @param serie the time series 
	 @return  forecastability level 
	*/
	public final double HowForecastable(ArrayList<Double> serie)
	{
		return HowForecastable(serie, true);
	}

	/**  Level of forecastability 
	 @param serie the time series 
	 @param trunk if it should be trunked 
	 @return  forecastability level 
	*/
	public final double HowForecastable(ArrayList<Double> serie, boolean trunk)
	{
		if (serie.size() < fcstHorizon)
		{
			throw new RuntimeException("El_horizonte_de_prediccion_no_puede_ser_mayor_que_el_largo_de_la_serie");
		}

		double[] real = new double[fcstHorizon];
		int size = serie.size();
		double maxError = 0.0;
		for (int i = 0;i < fcstHorizon;i++)
		{
			real[i] = (double)serie.get(size - fcstHorizon + i);
			maxError += real[i];
		}
		ArrayList<Double> datos = new ArrayList<Double>();
		datos.addAll(serie);
		datos.removeRange(size - fcstHorizon, fcstHorizon + size - fcstHorizon);
		tdnn = new TDNN(mw);
		tdnn.LoadSerie(datos, 0);
		tdnn.Train(ciclos);
		double[] fcst = tdnn.Process(ciclos,fcstHorizon);
		//Console.Write("\nData: "); foreach(double val in datos) { Console.Write(val.ToString("0.00") + "\t"); }
		//Console.Write("\nReal: "); foreach(double val in real) { Console.Write(val.ToString("0.00")+ "\t"); }
		//Console.Write("\nFcst: "); foreach(double val in fcst) { Console.Write(val.ToString("0.00")+ "\t"); }
		double error = GetLagError(fcst, real, 1);
		if (!trunk)
		{
			return error;
		}
		double relError = error / maxError;
		if (relError > 1)
		{
			relError = 1;
		}
		double forecastable = 1.0 - relError;
		//Console.WriteLine("\nError : " + error.ToString("0.00")+ " - RelErr: " + relError.ToString("0.00%"));
		//Console.WriteLine("Fcstbl: " + forecastable.ToString("0.00%"));
		return forecastable;
	}

	//endregion

	//region Private Methods

	private boolean HasUltPerZero(ArrayList<Double> serie)
	{
		for (int i = serie.size() - nUltPeriodos; i < serie.size();i++)
		{
			if ((double)serie.get(i) != 0)
			{
				return false;
			}
		}
		return true;
	}

	private double GetFc(ArrayList<Double> serie)
	{
		int totNoZero = 0;
		for (double val : serie)
		{
			if (val != 0)
			{
				totNoZero++;
			}
		}
		return (double)totNoZero / (double)serie.size();
	}

	private double GetLagError(double[] fcst, double[] real, int lag)
	{
		double errNorm = GetError(fcst, real, lag, lag, fcst.length - lag);
		double errLag = GetError(fcst, real, 0, lag, fcst.length - lag);
		//Console.WriteLine("\nErrNor : " + errNorm.ToString("0.00"));
		//Console.WriteLine("\nErrLag : " + errLag.ToString("0.00"));
		return Math.min(errNorm,errLag);
	}

	private double GetError(double[] fcst, double[] real)
	{
		return GetError(fcst, real, 0, 0, fcst.length);
	}

	private double GetError(double[] fcst, double[] real, int index1, int index2, int nDatos)
	{
		double error = 0.0;
		for (int i = 0;i < nDatos;i++)
		{
			error += Math.abs(fcst[index1 + i] - real[index2 + i]);
		}
		return error;
	}

	//endregion

	//region Metodos Debugging

	/**  Process from an input file, and generate an output file 
	 @param inPath input file path 
	 @param outPath output file path 
	*/
	public final void Process(String inPath, String outPath)
	{
		ArrayList<ArrayList<Double>> series = new ArrayList<ArrayList<Double>>();
		ArrayList<Double> real = new ArrayList<Double>();
		ArrayList<Double> fcst = new ArrayList<Double>();
		java.io.InputStreamReader sr = new java.io.InputStreamReader(inPath);
		int subPred = 0;
		int sobPred = 0;
		int okPred = 0;
		String[] tokens;
		char separador = '\t';
		ArrayList<Double> serie;
		double val;
		String line = sr.ReadLine();
		while (line != null)
		{
			serie = new ArrayList<Double>();
			tokens = line.split(java.util.regex.Pattern.quote(separador.toString()), -1);
			//precio = Convert.ToDouble(tokens[0]); 
			for (int i = 1;i < 37;i++)
			{
				val = Double.parseDouble(tokens[i]);
				serie.add(val);
			}
			series.add(serie);
			real.add(Short.parseShort(tokens[37]));
			line = sr.ReadLine();
		}
		sr.close();


		int fcstable;
		for (ArrayList<Double> ser : series)
		{
			fcstable = (IsForecastable(ser)? 1 : 0);
			fcst.add(fcstable);
		}

		java.io.OutputStreamWriter sw = new java.io.OutputStreamWriter(outPath);
		sw.write("RESULTADOS PREDECIBILIDAD" + System.lineSeparator());
		sw.write("" + System.lineSeparator());
		sw.write("FT\tRL\tSERIE" + System.lineSeparator());
		sw.write("" + System.lineSeparator());

		for (int i = 0;i < series.size();i++)
		{
			if (fcst.get(i).toString().equals("1") && real.get(i).toString().equals("0"))
			{
				sobPred++;
			}
			else if (fcst.get(i).toString().equals("0") && real.get(i).toString().equals("1"))
			{
				subPred++;
			}
			else
			{
				okPred++;
			}

			sw.write(fcst.get(i) + "\t" + real.get(i) + "\t");
			serie = (ArrayList<Double>)series.get(i);
			for (double valor : serie)
			{
				sw.write(valor + "\t");
			}
			sw.write("" + System.lineSeparator());
		}
		sw.write("" + System.lineSeparator());
		sw.write("RESULTADOS: " + System.lineSeparator());
		sw.write("" + System.lineSeparator());
		sw.write("Prediccion OK  : " + okPred + System.lineSeparator());
		sw.write("SubPrediccion  : " + subPred + System.lineSeparator());
		sw.write("SobrePrediccion: " + sobPred + System.lineSeparator());
		sw.write("Total          : " + series.size() + System.lineSeparator());
		sw.close();
	}

	//endregion
}