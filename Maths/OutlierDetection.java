package NumCalc;

import ClassicalStat.*;
import java.util.*;

/**  Class for outlier detection 
*/
public class OutlierDetection
{

	//region Fields

	private double[] xArr = {90.0, 90.5, 91.0, 91.5, 92.0, 92.5, 93.0, 93.5, 94.0, 94.5, 95.0, 95.5, 96.0, 96.5, 97.0, 97.5, 98.0, 98.5, 99.0, 99.3, 99.5, 99.7, 99.9};
	private double[] yArr = {1.88599968, 1.95, 2.05, 2.15, 2.24, 2.3, 2.41, 2.4, 2.55, 2.65, 2.919998884, 3.1, 3.3, 3.5, 3.8, 4.302995682, 4.8, 5.6, 6.964999676, 8.2, 9.925000668, 13.0, 22.32700014};

	private Splines sp;
	private Statistics stat;

	//endregion

	//region Constructor

	/**  Constructor 
	*/
	public OutlierDetection()
	{
		this.sp = new Splines();
		this.stat = new Statistics();
	}

	//endregion

	//region Grubbs Test

	/**  Test que indica si hay al menos un outlier en la serie dada 
	 @param datos data 
	 @param prob probability 
	 @return  if there are outlier 
	*/
	public final boolean GrubbsTest(List<Double> datos, double prob)
	{
		int minLenghtForOutliers = 3;
		if (datos.size() <= minLenghtForOutliers)
		{
			return false;
		}
		//determinacion del valor G = máximo valor (val-m)/stDev
		double G = 0.0;
		double mean = stat.Mean(datos);
		double stDev = stat.StDev(datos);
		double g;
		for (double val : datos)
		{
			g = (val - mean) / stDev;
			if (g > G)
			{
				G = g;
			}
		}
		double k = GetCriticalValue(prob);

		//determinacion del valor critico (cv);
		int n = datos.size();
		double cv = ((n - 1) / Math.sqrt(n)) * Math.sqrt(Math.pow(k, 2) / (n - 2 + Math.pow(k, 2)));

		//comparacion final
		//Console.WriteLine("G = "+ G + " , cv = " + cv + " , n = " + n + " , prob = " + prob);
		return G > cv;
	}

	/**  Test de identificacion de outliers individual
	 @param datos data 
	 @param index index 
	 @param prob probability 
	 @param margen grubbs margin 
	 @return  if it is an outlier 
	*/
	public final boolean GrubbsTest(List<Double> datos, int index, double prob, double margen)
	{

		int minLenghtForOutliers = 3;
		if (datos.size() <= minLenghtForOutliers)
		{
			return false;
		}

		//determinacion del valor G = (val-m)/stDev
		double mean = stat.Mean(datos);
		double stDev = stat.StDev(datos);
		double G = (datos.get(index) - mean) / stDev;

		double k = GetCriticalValue(prob);
		//Console.WriteLine("df = " + df);

		//determinacion del valor critico (cv);
		int n = datos.size();
		double cv = ((n - 1) / Math.sqrt(n)) * Math.sqrt(Math.pow(k, 2) / (n - 2 + Math.pow(k, 2)));

		//comparacion final
		return G / cv > 1 - margen / 100.0;
	}


	/**  Test de identificacion de outliers individual sobre Frecuencias 
	 @param frec frequency 
	 @param value value 
	 @param prob probability 
	 @param margen grubbs margin 
	 @return  if it is an outlier 
	*/
	public final boolean GrubbsTest(Frequencies frec, double value, double prob, double margen)
	{

		if (value < 0.0)
		{
			return true;
		}
		int minLenghtForOutliers = 3;
		if (frec.NValues <= minLenghtForOutliers)
		{
			return false;
		}

		//determinacion del valor G = (val-m)/stDev

		int nValues = frec.NValues + 1;
		double mean = frec.Mean;
		double stDev = frec.StDev;

		double G = (value - mean) / stDev;

		double k = GetCriticalValue(prob);
		//Console.WriteLine("df = " + df);

		//determinacion del valor critico (cv);
		double cv = ((nValues - 1) / Math.sqrt(nValues)) * Math.sqrt(Math.pow(k, 2) / (nValues - 2 + Math.pow(k, 2)));

		//comparacion final
		return G / cv > 1 - margen / 100.0;
	}

	/**  Cargar splines 
	*/
	public final void LoadSplines()
	{
		//determinacion del df para t de student n-2 grados de libertad
		ArrayList<Double> Xdata = new ArrayList<Double>();
		tangible.DoubleLists.addPrimitiveArrayToList(xArr, Xdata);
		ArrayList<Double> Ydata = new ArrayList<Double>();
		tangible.DoubleLists.addPrimitiveArrayToList(yArr, Ydata);
		sp = new Splines(Xdata, Ydata, true);
	}

	private double GetCriticalValue(double prob)
	{

		prob = Math.round(prob * Math.pow(10, 1)) / Math.pow(10, 1);

		double x = prob;
		double y = 0.0;
		double x1 = 0;
		double y1 = 0;
		double x2 = 0;
		double y2 = 0;
		for (int j = 0;j < xArr.length;j++)
		{
			if (xArr[j] == x)
			{
				y = yArr[j];
				return y;
			}
			else if (xArr[j] > x)
			{
				x1 = xArr[j - 1];
				x2 = xArr[j];
				y1 = yArr[j - 1];
				y2 = yArr[j];
				if (sp == null)
				{
					LoadSplines();
				}
				y = sp.Interpolar(x, x1, y1, x2, y2);
				return y;
			}
		}
		throw new RuntimeException(String.format(Strings.No_pude_interpolarse_el_valor_0, prob));
	}

	//endregion

	//region Studentized deletion residuals

	//region Public Methods

	/**  Standardized residuals 
	 @param serie list of data values of the time series 
	 @param model list of data values of the model 
	 @return  list of standardized residuals 
	*/
	public final ArrayList<Double> StdRes(List<Double> serie, List<Double> model)
	{
		if (serie.size() != model.size())
		{
			throw new RuntimeException(Strings.Serie_and_model_must_have_the_same_size);
		}

		double s = Math.sqrt(Mse(serie, model));
		ArrayList<Double> stdres = new ArrayList<Double>();
		for (int i = 0;i < serie.size();i++)
		{
			stdres.add((serie.get(i) - model.get(i)) / s);
		}
		return stdres;
	}

	/**  Studentized residuals 
	 @param serie list of data values of the time series 
	 @param model list of data values of the model 
	 @return  list of studentized residuals 
	*/
	public final ArrayList<Double> StudRes(List<Double> serie, List<Double> model)
	{
		if (serie.size() != model.size())
		{
			throw new RuntimeException(Strings.Serie_and_model_must_have_the_same_size);
		}

		ArrayList<Double> studres = new ArrayList<Double>();
		ArrayList<Double> stdres = StdRes(serie, model);
		ArrayList<Double> hat = Hat(serie, model);
		for (int i = 0;i < serie.size();i++)
		{
			studres.add(stdres.get(i) / Math.sqrt(1 - hat.get(i)));
		}
		return studres;
	}

	/**  Studentized deletion residuals 
	 @param serie list of data values of the time series 
	 @param model list of data values of the model 
	 @return  list of studentized deletion residuals 
	*/
	public final ArrayList<Double> StudDelRes(List<Double> serie, List<Double> model)
	{
		if (serie.size() != model.size())
		{
			throw new RuntimeException(Strings.Serie_and_model_must_have_the_same_size);
		}

		ArrayList<Double> studDelRes = new ArrayList<Double>();
		ArrayList<Double> r = StudRes(serie, model);
		ArrayList<Double> hat = Hat(serie, model);
		double n = serie.size();
		for (int i = 0;i < n;i++)
		{
			if (n - 2 - Math.pow(r.get(i), 2) < 0)
			{
				studDelRes.add(Double.MAX_VALUE);
			} //TODO: revisar
			else
			{
				studDelRes.add(r.get(i) * Math.sqrt((n - 3) / (n - 2 - Math.pow(r.get(i), 2))));
			}
		}
		//for(int i=0;i<n;i++) { studDelRes.Add((serie[i]-model[i]) / (1-hat[i])); } //TODO: revisar si el numerador es r o el raw residual
		//for(int i=0;i<n;i++) { studDelRes.Add(r[i]/ (1-hat[i])); }
		return studDelRes;
	}

	/**  Limit for outlier filtering with studentized deletion residuals 
	 @param serie list of data values of the time series 
	 @param model list of data values of the model 
	 @param alpha significance level 
	 @return  the limit value 
	*/
	public final double StudDelResLimit(List<Double> serie, List<Double> model, double alpha)
	{
		return Stat.TStudent_quantil(alpha, serie.size(), true);
	}

	/**  List of indexes of detected outliers 
	 @param serie list of data values of the time series 
	 @param model list of data values of the model 
	 @param alpha significance level 
	 @return  list of indexes of outliers 
	*/
	public final ArrayList<Integer> StudDelResOutliers(List<Double> serie, List<Double> model, double alpha)
	{
		ArrayList<Integer> outliers = new ArrayList<Integer>();
		ArrayList<Double> studDelRes = StudDelRes(serie, model);
		double limit = StudDelResLimit(serie, model, alpha);
		for (int i = 0;i < studDelRes.size();i++)
		{
			if (Math.abs(studDelRes.get(i)) > limit)
			{
				outliers.add(i);
			}
		}
		return outliers;
	}

	/**  If it is an outlier according to standardized deletion resiuduals 
	 @param serie time series 
	 @param model model used 
	 @param alpha significance 
	 @param value value to test 
	 @param modelValue model value for test 
	 @return 
	*/
	public final boolean IsStudDelResOutlier(List<Double> serie, List<Double> model, double alpha, double value, double modelValue)
	{
		serie.add(value);
		model.add(modelValue);
		ArrayList<Double> studDelRes = StudDelRes(serie, model);
		double limit = StudDelResLimit(serie, model, alpha);
		if (Math.abs(studDelRes.get(studDelRes.size() - 1)) > limit)
		{
			return true;
		}
		return false;
	}

	public final boolean IsStudDelResOutlier(List<Double> serie, List<Double> model, ArrayList<Double> studDelRes, double alpha, int index)
	{
		double limit = StudDelResLimit(serie, model, alpha);
		if (Math.abs(studDelRes.get(index)) > limit)
		{
			return true;
		}
		return false;
	}

	/**  Leverage values 
	 @param serie list of data values of the time series 
	 @param model list of data values of the model 
	 @return  list of leverage values 
	*/
	public final ArrayList<Double> Leverage(List<Double> serie, List<Double> model)
	{
		if (serie.size() != model.size())
		{
			throw new RuntimeException(Strings.Serie_and_model_must_have_the_same_size);
		}
		ArrayList<Double> hat = Hat(serie, model);
		ArrayList<Double> lev = new ArrayList<Double>();
		for (int i = 0;i < serie.size();i++)
		{
			lev.add(1 - hat.get(i));
		}
		return lev;
	}

	//endregion 

	//region Private Methods

	private double Sse(List<Double> serie, List<Double> model)
	{
		if (serie.size() != model.size())
		{
			throw new RuntimeException(Strings.Serie_and_model_must_have_the_same_size);
		}

		double sse = 0.0;
		for (int i = 0;i < serie.size();i++)
		{
			sse = sse + Math.pow(serie.get(i) - model.get(i), 2);
		}
		return sse;
	}

	private double Mse(List<Double> serie, List<Double> model)
	{
		if (serie.size() != model.size())
		{
			throw new RuntimeException(Strings.Serie_and_model_must_have_the_same_size);
		}
		double sse = Sse(serie, model);
		return sse / (serie.size() - 2);
	}

	private ArrayList<Double> Hat(List<Double> serie, List<Double> model)
	{
		if (serie.size() != model.size())
		{
			throw new RuntimeException(Strings.Serie_and_model_must_have_the_same_size);
		}

		double sse = Sse(serie, model);
		ArrayList<Double> hat = new ArrayList<Double>();
		for (int i = 0;i < serie.size();i++)
		{
			hat.add(1 / serie.size() + Math.pow(serie.get(i) - model.get(i),2) / sse);
		}
		return hat;
	}

	//endregion

	//endregion

}