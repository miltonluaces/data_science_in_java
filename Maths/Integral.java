package NumCalc;

import java.util.*;

//region Imports


//endregion


/**  Class for simpson numerical integration 
*/
public class Integral
{

	//region Fields

	private Hashtable data;
	private ArrayList<Double> coefcs;
	private Splines sp;
	private TipoRegresion tipo = TipoRegresion.values()[0];
	private Polynomic poly;

	//endregion

	//region Constructors

	/**  Constructor por defecto. Asigna Regresion de Cubic Splines 
	*/
	public Integral()
	{
		tipo = TipoRegresion.LinearSpline;
		poly = new Polynomic();
	}

	/**  Constructor con parametro para asignar tipo de regresion 
	 @param tipo type of regression 
	*/
	public Integral(TipoRegresion tipo)
	{
		this.tipo = tipo;
	}

	//endregion

	//region Public Methods

	/**  Carga de Datos desde ArrayList  
	 @param x independent variable 
	 @param y dependent variable (data) 
	*/
	public final void LoadData(ArrayList<Double> x, ArrayList<Double> y)
	{
		switch (tipo)
		{
			case Aritmetica:
				break;
			case Polinomica:
				coefcs = poly.Regression(x, y, x.size() - 1);
				break;
			case LinearSpline:
				sp = new Splines(x, y, true);
				break;
			case CubicSpline:
				sp = new Splines(x, y, false);
				break;
		}
	}

	/**  Carga de Datos desde Hash 
	 @param data the hashtable 
	*/
	public final void LoadData(Hashtable data)
	{
		this.data = data;
	}

	/**  Obtain the splines 
	 @param min min value 
	 @param max max value 
	 @param step step for discretization 
	 @return  list of splines 
	*/
	public final ArrayList<Double> GetSplines(int min, int max, int step)
	{
		ArrayList<Double> splines = new ArrayList<Double>();
		for (int i = min;i <= max;i = i + step)
		{
			splines.add(fSpline(i));
		}
		return splines;
	}

	/**  Integral con 10 iteraciones. Por defecto 
	 @param a left integration limit 
	 @param b right integration limit 
	 @return  value of calculated integral 
	*/
	public final double IntegralDef(double a, double b)
	{
		return IntegralDef(a, b, 100);
	}

	/**  
	 Implementaci�n del M�todo de Simpson
	 Parte el intervalo en 2 * h y luego pondera y1 + 4y2 + y3 
	 int = h/3 (y1 + 4y2 + y3)
	 
	 @param a limite de integracion izq 
	 @param b limite de integracion derecho 
	 @param n numero de columnas 
	 @return  value of calculated integration 
	*/
	public final double IntegralDef(double a, double b, int n)
	{
		//conversion a par
		if (n % 2 == 1)
		{
			n++;
		}

		double h = (b - a) / n;

		//extremos
		double suma = f(a) + f(b);
		//impares
		for (int i = 1;i < n;i += 2)
		{
			suma += 4 * f(a + i * h);
		}
		//pares
		for (int i = 2;i < n;i += 2)
		{
			suma += 2 * f(a + i * h);
		}
		return suma * h / 3;
	}

	//endregion

	//region Private Methods

	/**  Default interpolation 
	 @param x value to interpolate (independent variable) 
	 @return  interpolated value 
	*/
	public final double f(double x)
	{
		switch (tipo)
		{
			case Aritmetica:
				return fAritm(x);
			case Polinomica:
				return fRegPol(x);
			case LinearSpline:
				return fSpline(x);
			case CubicSpline:
				return fSpline(x);
		}
		return -Double.MAX_VALUE;
	}

	/**  Calculate n inverse values for y (not necesarily biyective function) 
	 @param y functional value 
	 @return  n x inverse values 
	*/
	public final ArrayList<Double> Inversas(double y)
	{
		return sp.Inversas(y);
	}

	/**  Calculate positive intervals 
	 @param y functional value 
	 @return  list of points that define the intervals 
	*/
	public final ArrayList<Splines.Point> IntervalosPositivos(double y)
	{
		return sp.IntervalosPositivos(y);
	}


	/**  Spline Interpolation 
	 @param x value to interpolate (independent variable) 
	 @return  interpolation 
	*/
	public final double fSpline(double x)
	{
		return sp.Interpolar(x);
	}

	/**  Polynomic regression interpolation 
	 @param x value to interpolate (independent variable) 
	 @return  interpolation 
	*/
	public final double fRegPol(double x)
	{
		double res = 0.0;
		for (int g = 0;g < coefcs.size();g++)
		{
			res += (double)coefcs.get(g) * Math.pow(x, g);
		}
		return res;
	}

	/**  Arithmetic interpolation 
	 @param x value to interpolate (independent variable) 
	 @return  interpolation 
	*/
	private double fAritm(double x)
	{

		if (data.get(x) != null)
		{
			return (double)data.get(x);
		}
		else
		{
			double difDef = Double.MAX_VALUE;
			double difExc = Double.MAX_VALUE;
			double xDef = Double.MAX_VALUE;
			double xExc = Double.MAX_VALUE;
			for (double xAct : data.keySet())
			{
				if (xAct < x && Math.abs(xAct - x) < difDef)
				{
					difDef = Math.abs(xAct - x);
					xDef = xAct;
				}
				else if (xAct > x && Math.abs(xAct - x) < difExc)
				{
					difExc = Math.abs(xAct - x);
					xExc = xAct;
				}
			}
			return ((double)data.get(xExc) + (double)data.get(xExc)) / 2.00;
		}
	}

	//endregion

	//region Enum TipoRegresion

	/**  Regression type 
	*/
	public enum TipoRegresion
	{
		/**  Arithmetic regression 
		*/
		Aritmetica,
		/**  Polynomic regression 
		*/
		Polinomica,
		/**  Linear spline 
		*/
		LinearSpline,
		/**  Cubic b-splines 
		*/
		CubicSpline;

		public int getValue()
		{
			return this.ordinal();
		}

		public static TipoRegresion forValue(int value)
		{
			return values()[value];
		}
	}

	//endregion
}