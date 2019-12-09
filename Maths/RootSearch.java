package NumCalc;

//region Imports

import ClassicalStat.*;
import java.util.*;

//endregion


/**  Repository of algorithms for function roots search 
*/
public class RootSearch  {

	//region Fields

	private double epsilon;
	private double epsilonInterval;
	private int iterations;
	private int maxIterations;

	private HashMap<String, Double> fCol;
	private HashMap<String, Double> fdxCol = null;
	private HashMap<String, Double> fdx2Col = null;

	private RandomGen randGen;

	//endregion

	//region Public Functions

	/**  Function 
	 @param x independent variable 
	 @return  dependent variable 
	*/
	public final double f(double x)
	{
		return fCol.get((new Double(x)).toString("0.00"));
	}

	/**  First derivative Function 
	 @param x independent variable 
	 @return  dependent variable 
	*/
	public final double fdx(double x)
	{
		return fdxCol.get((new Double(x)).toString("0.00"));
	}

	/**  Second derivative Function 
	 @param x independent variable 
	 @return  dependent variable 
	*/
	public final double fdx2(double x)
	{
		return fdx2Col.get((new Double(x)).toString("0.00"));
	}

	//endregion

	//region Public Delegates

	/** 
	 <p>Función de la cual se desea encontrar el cero o ceros</p>
	 <p>Es un delegado, por tanto con el delegado encontramos el valor de la función en x</p>
	 
	 @param x Valor para el cual deseamos obtener la funion
	*/
	@FunctionalInterface
	public interface function
	{
		double invoke(double x);
	}
	/** 
	 <p>Derivada primera de la función de la cual se desea encontrar el cero o ceros</p>
	 <p>Es un delegado, por tanto con el delegado encontramos el valor de la derivada primera en x</p>
	 
	 @param x Valor para el cual deseamos obtener la derivada primera de la funcion
	*/
	@FunctionalInterface
	public interface derivative
	{
		double invoke(double x);
	}
	/** 
	 <p>Derivada segunda de la función de la cual se desea encontrar el cero o ceros</p>
	 <p>Es un delegado, por tanto con el delegado encontramos el valor de la derivada segunda en x</p>
	 
	 @param x Valor para el cual deseamos obtener la derivada segunda de la funcion
	*/
	@FunctionalInterface
	public interface derivative2
	{
		double invoke(double x);
	}

	//endregion

	//region Constructor

	/** 
	 <p>Constructor de la clase</p>
	 
	 @param epsilon Error maximo de la raiz o cero buscado
	 @param epsilonInterval Porcentaje del rango de intervalo sobre el valor
	 @param maxIterations Maximo de iteraciones
	*/
	public RootSearch(double epsilon, double epsilonInterval, int maxIterations)
	{
		this.epsilon = epsilon;
		this.epsilonInterval = epsilonInterval;
		this.maxIterations = maxIterations;
		this.randGen = new RandomGen();
	}

	//endregion

	//region Properties

	/**  number of iterations allowed 
	*/
	public final int getIterations()
	{
		return iterations;
	}

	/**  Epsilon allowed for solution validation 
	*/
	public final double getEpsilon()
	{
		return epsilon;
	}
	public final void setEpsilon(double value)
	{
		epsilon = value;
	}

	//endregion

	//region Range Methods

	public final HashMap<Double, Double> RangeSearch(function f, double xMin, double xMax, double yMin, double yMax, double step, double eps, tangible.RefObject<Integer> its, int maxIt)
	{
		int hits = 0;
		int it = 0;
		double relEpsX = 0.05;
		int nRanges = -1;
		double x, y, xIni, xEnd, yIni, yEnd, keY;
		randGen.Reset();
		HashMap<Double, Double> yxValues = new HashMap<Double, Double>();
		try
		{
			//initial values
			yxValues.put(yMin, xMin);
			y = yMin + step;
			while (y < yMax)
			{
				yxValues.put(y, -1);
				y = Math.round((y + step) * Math.pow(10, 1)) / Math.pow(10, 1);
			}
			yxValues.put(yMax, xMax);

			tangible.RefObject<Integer> tempRef_it = new tangible.RefObject<Integer>(it);
			hits += FillValues((double x) -> f(x), yxValues, xMin, xMax, tempRef_it, maxIt, 1);
		it = tempRef_it.argValue;
			its.argValue += it;
			xIni = xMin;
			yIni = yMin;
			int nMiss = 0;
			ArrayList<Range> ranges = new ArrayList<Range>();
			int wIts = 0;
			int maxWIts = 10;
			while (true)
			{
				for (double key : yxValues.keySet())
				{
					y = key;
					x = yxValues.get(key);
					if (x == -1)
					{
						nMiss++;
					}
					else
					{
						if (nMiss == 0)
						{
							xIni = x;
							yIni = y;
						}
						else
						{
							xEnd = x;
							yEnd = y;
							ranges.add(new Range(xIni, xEnd, yIni, yEnd, nMiss));
							nMiss = 0;
							xIni = x;
							yIni = y;
						}
					}
				}
				if (ranges.isEmpty())
				{
					break;
				}
				for (Range range : ranges)
				{
					if (range.nValues == 1)
					{
						keY = Math.round((range.yIni + step) * Math.pow(10, 1)) / Math.pow(10, 1);
						if (keY <= yMax)
						{
							yxValues.put(keY, (range.xIni + range.xEnd) / 2.0);
						}
						hits++;
					}
					else
					{
						if (wIts > maxWIts || ranges.size() == nRanges || (range.xEnd - range.xIni) / range.xEnd < relEpsX)
						{
							y = range.yIni;
							x = (range.xIni + range.xEnd) / 2.0;
							while (y <= range.yEnd)
							{
								y = y + step;
								keY = Math.round(y * Math.pow(10, 1)) / Math.pow(10, 1);
								if (keY == 100.0 && !yxValues.containsKey(100.0))
								{
									yxValues.put(keY, yxValues.get(99.9));
									hits++;
								}
								if (yxValues.containsKey(keY) && yxValues.get(keY).equals(-1))
								{
									yxValues.put(keY, x);
									hits++;
								}
							}
						}
						else
						{
							it = 0;
							tangible.RefObject<Integer> tempRef_it2 = new tangible.RefObject<Integer>(it);
							hits += FillValues((double x) -> f(x), yxValues, range.xIni, range.xEnd, tempRef_it2, range.nValues * 10, 1);
						it = tempRef_it2.argValue;
							its.argValue += it;

						}
					}
				}
				if (yxValues.get(99.9).equals(-1))
				{
					yxValues.put(99.9, yxValues.get(100));
					hits++;
				}
				nRanges = ranges.size();
				ranges.clear();

				wIts++;
			}
		}
		catch (RuntimeException ex)
		{
			throw new RuntimeException("Error in range calculation: " + ex.getMessage(), ex);
		}
		double prevValue = -1;
		HashMap<Double, Double> monotoneYxValues = new HashMap<Double, Double>();
		for (double key : yxValues.keySet())
		{
		   if (yxValues.get(key).compareTo(prevValue) >= 0)
		   {
			   monotoneYxValues.put(key, yxValues.get(key));
		   }
		   else
		   {
			   monotoneYxValues.put(key, prevValue);
		   }
		   prevValue = yxValues.get(key);
		}
		return monotoneYxValues;
	}

	private int FillValues(function f, HashMap<Double, Double> yxValues, double xMin, double xMax, tangible.RefObject<Integer> it, int maxIt, int dec)
	{
		if (xMin == xMax)
		{
			return 0;
		}
		double x, y;
		int hits = 0;
		for (it.argValue = 0; it.argValue < maxIt; it.argValue++)
		{
			x = randGen.NextDouble(xMin, xMax);
			y = Math.round((f.invoke(x) * 100) * Math.pow(10, dec)) / Math.pow(10, dec);
			if (yxValues.containsKey(y) && yxValues.get(y).equals(-1))
			{
				yxValues.put(y, x);
				hits++;
			}
		}
		return hits;
	}


	//endregion

	//region Bisection Methods

	/**  Bisection root search numerical method 
	 @param f function 
	 @param d first derivative 
	 @param min min value of the interval 
	 @param max max value of the interval 
	 @param val value to add (0 if pure root search) 
	 @param eps epsilon for solution validation 
	 @param maxIt maximum of iterations 
	 @return  value of independent variable 
	*/
	public final double MonotoneBisection(function f, derivative d, double min, double max, double val, double eps, tangible.RefObject<Integer> it, int maxIt)
	{
		double dMin = d.invoke(min);
		double dMax = d.invoke(max);
		if (dMin * dMax < 0)
		{
			throw new IllegalArgumentException(Strings.This_method_is_only_for_monotone_intervals);
		}
		return MonotoneBisectionRec((double x) -> f(x), min, max, val, (dMin > 0), eps, it, maxIt);
	}

	/**  Bisection root search numerical method 
	 @param f function 
	 @param ascendant if the function is monotone, ascendant or not 
	 @param min min value of the interval 
	 @param max max value of the interval 
	 @param val value to add (0 if pure root search) 
	 @param eps epsilon for solution validation 
	 @param maxIt maximum of iterations 
	 @return  value of independent variable 
	*/
	public final double MonotoneBisection(function f, boolean ascendant, double min, double max, double val, double eps, tangible.RefObject<Integer> it, int maxIt)
	{
		return MonotoneBisectionRec((double x) -> f(x), min, max, val, ascendant, eps, it, maxIt);
	}

	/**  Bisection root search numerical method 
	 @param f function 
	 @param min min value of the interval 
	 @param max max value of the interval 
	 @param val value to add (0 if pure root search) 
	 @param ascendant if the function is monotone, ascendant or not 
	 @param eps epsilon for solution validation 
	 @param maxIt maximum of iterations 
	 @param it current iteration 
	 @return  value of independent variable 
	*/
	private double MonotoneBisectionRec(function f, double min, double max, double val, boolean ascendant, double eps, tangible.RefObject<Integer> it, int maxIt)
	{
		it.argValue++;
		double inter = min + (max - min) / 2.0;
		if (it.argValue > maxIt || inter == min || inter == max)
		{
			//Console.WriteLine("Error.More than 30 iterations.");
			return inter;
		}
		double fInter = f.invoke(inter);
		if (Math.abs(val - fInter) < eps)
		{
			return inter;
		}
		if (ascendant)
		{
		   if (val < fInter)
		   {
			   return MonotoneBisectionRec((double x) -> f(x), min, inter, val, ascendant, eps, it, maxIt);
		   }
		   if (val > fInter)
		   {
			   return MonotoneBisectionRec((double x) -> f(x), inter, max, val, ascendant, eps, it, maxIt);
		   }
		}
		else
		{
			if (val > fInter)
			{
				return MonotoneBisectionRec((double x) -> f(x), min, inter, val, ascendant, eps, it, maxIt);
			}
			if (val < fInter)
			{
				return MonotoneBisectionRec((double x) -> f(x), inter, max, val, ascendant, eps, it, maxIt);
			}
		}
		//Console.WriteLine("Error inter = " + inter);
		return inter;
	}

	//endregion

	//region Newton Raphson Methods

	/** 
	 <p>Metodo Newton Raphson para encontrar una raiz o cero de una función</p>
	 <p>Si <paramref name="val"/> es cero entonces encuantra la raiz de la funcion <paramref name="f"/>,</p>
	 <p>si no, encuentra la inversa de la función <paramref name="f"/> en <paramref name="val"/></p>
	 
	 @param f funcion mediante la cual se calcula el valor de dicha funcion
	 @param d funcion mediante la cual se calcula el valor de la derivada primera de dicha funcion
	 @param d2 funcion mediante la cual se calcula el valor de la derivada segunda de dicha funcion
	 @param min Valor minimo del intervalo en el cual busca la raiz o cero
	 @param max Valor maximo del intervalo en el cual busca la raiz o cero
	 @param val Valor para el cual se desea buscar la inversa de la función <paramref name="f"/>
	 @return  root value 
	*/
	public final double NewtonRaphsonOneRoot(function f, derivative d, derivative2 d2, double min, double max, double val)
	{
		if ((f.invoke(min) - val) * f.invoke(max) - val > 0)
		{
			throw new IllegalArgumentException(Strings.No_zero_between_these_values);
		}
		double nr = 0;
		double x;
		int i = 0;
		while (nr == 0)
		{
			x = GetFourierValue((double x) -> f(x), d2, min, max, val);
			nr = NewtonRaphson((double x) -> f(x), d, x, val);
			if (i > maxIterations)
			{
				throw new RuntimeException(Strings.Exceed_the_maximun_iterations);
			}
			i++;
		}
		return nr;
	}

	/** 
	 <p>Metodo Newton Raphson para encontrar una raiz o cero de una función</p>
	 <p>Si <paramref name="val"/> es cero entonces encuantra la raiz de la funcion <paramref name="f"/>,</p>
	 <p>si no, encuentra la inversa de la función <paramref name="f"/> en <paramref name="val"/></p>
	 
	 @param f funcion mediante la cual se calcula el valor de dicha funcion
	 @param d funcion mediante la cual se calcula el valor de la derivada primera de dicha funcion
	 @param x independent variable 
	 @param val Valor para el cual se desea buscar la inversa de la función <paramref name="f"/>
	 @return  root value 
	*/
	public final double NewtonRaphson(function f, derivative d, double x, double val)
	{
		int i = 0;
		double gx = f.invoke(x) - val;
		while (Math.abs(gx) > epsilon)
		{
			if (d.invoke(x) == 0)
			{
				return 0;
			}
			x = x - gx / d.invoke(x);
			gx = f.invoke(x) - val;
			i++;
			if (i > maxIterations)
			{
				throw new RuntimeException(Strings.Exceed_the_maximun_iterations);
			}
		}
		return x;
	}

	//endregion

	//region Steffensen Acceleration - Aitken Iteration

	private double GetAitkenIteration(double xiMin2, double xiMin1, double xi)
	{
		double num = xi * xiMin2 - quad(xiMin1);
		double den = (xi - 2 * xiMin1 + xiMin2);
		double res = num / den;
		return (xi * xiMin2 - quad(xiMin1)) / (xi - 2 * xiMin1 + xiMin2);
	}

	/** 
	 <p>Metodo Steffenson, que es una forma mas rapida de encontrar una raiz o cero de una función</p>
	 <p>Si <paramref name="val"/> es cero entonces encuantra la raiz de la funcion <paramref name="f"/>,</p>
	 <p>si no, encuentra la inversa de la función <paramref name="f"/> en <paramref name="val"/></p>
	 
	 @param f funcion mediante la cual se calcula el valor de dicha funcion
	 @param d funcion mediante la cual se calcula el valor de la derivada primera de dicha funcion
	 @param x independent variable 
	 @param val Valor para el cual se desea buscar la inversa de la función <paramref name="f"/>
	 @return  root value 
	*/
	public final double SteffensenAcceleration(function f, derivative d, double x, double val)
	{
		double gxAnt = Double.MAX_VALUE;
		epsilon = 0.005;
		double gx = f.invoke(x) - val;
		int i = 0;
		double xi = -1;
		double xiMin1 = -1;
		double xiMin2 = -1;

		while (Math.abs(gx) > epsilon)
		{
			if (d.invoke(x) == 0)
			{
				return 0;
			}
			if (i < 4 || i % 2 == 0 || xi - 2 * xiMin1 + xiMin2 == 0)
			{
				x = x - gx / d.invoke(x);
			}
			else
			{
				x = GetAitkenIteration(xiMin2, xiMin1, xi);
			}
			gx = f.invoke(x) - val;
			if (Math.abs(gx) > Math.abs(gxAnt))
			{
				return xi;
			}
			gxAnt = gx;
			xiMin2 = xiMin1;
			xiMin1 = xi;
			xi = x;
			i++;
			if (i > maxIterations)
			{
				return 0;
			}
		}
		iterations = i;
		return x;
	}

	/** 
	 <p>Metodo Steffenson, que es una forma mas rapida de encontrar una raiz o cero de una función</p>
	 <p>Si <paramref name="val"/> es cero entonces encuantra la raiz de la funcion <paramref name="f"/>,</p>
	 <p>si no, encuentra la inversa de la función <paramref name="f"/> en <paramref name="val"/></p>
	 
	 @param f funcion mediante la cual se calcula el valor de dicha funcion
	 @param d funcion mediante la cual se calcula el valor de la derivada primera de dicha funcion
	 @param d2 funcion mediante la cual se calcula el valor de la derivada segunda de dicha funcion
	 @param min Valor minimo del intervalo en el cual busca la raiz o cero
	 @param max Valor maximo del intervalo en el cual busca la raiz o cero
	 @param val Valor para el cual se desea buscar la inversa de la función <paramref name="f"/>
	 @return  root value 
	*/
	public final double SteffensenAccOneRoot(function f, derivative d, derivative2 d2, double min, double max, double val)
	{
		if ((f.invoke(min) - val) * f.invoke(max) - val > 0)
		{
			throw new IllegalArgumentException(Strings.No_zero_between_these_values);
		}
		double st = 0;
		double x;
		int i = 0;
		while (st == 0)
		{
			x = GetFourierValue((double x) -> f(x), d2, min, max, val);
			st = SteffensenAcceleration((double x) -> f(x), d, x, val);
			i++;
			if (i > maxIterations)
			{
				throw new RuntimeException(Strings.Exceed_the_maximun_iterations);
			}
		}
		//Console.WriteLine("val = " + val + " , st = " + st + " , it = " + i);
		return st;
	}

	/** 
	 <p>Metodo Steffenson, que es una forma mas rapida de encontrar una raiz o cero de una función</p>
	 <p>Si <paramref name="val"/> es cero entonces encuantra la raiz de la funcion <paramref name="f"/>,</p>
	 <p>si no, encuentra la inversa de la función <paramref name="f"/> en <paramref name="val"/></p>
	 
	 @param f funcion mediante la cual se calcula el valor de dicha funcion
	 @param d funcion mediante la cual se calcula el valor de la derivada primera de dicha funcion
	 @param d2 funcion mediante la cual se calcula el valor de la derivada segunda de dicha funcion
	 @param min Valor minimo del intervalo en el cual busca la raiz o cero
	 @param max Valor maximo del intervalo en el cual busca la raiz o cero
	 @param val Valor para el cual se desea buscar la inversa de la función <paramref name="f"/>
	 @param maxIterations maximum of iterations 
	 @return  root value 
	*/
	public final double SteffensenAccOneRoot(function f, derivative d, derivative2 d2, double min, double max, double val, int maxIterations)
	{
		if ((f.invoke(min) - val) * f.invoke(max) - val > 0)
		{
			throw new RuntimeException(Strings.No_zero_between_these_values);
		}
		double st = 0;
		double x;
		int i = 0;
		while (st == 0)
		{
			x = GetFourierValue((double x) -> f(x), d2, min, max, val);
			st = SteffensenAcceleration((double x) -> f(x), d, x, val);
			i++;
			if (i > maxIterations)
			{
				throw new RuntimeException(Strings.Exceed_the_maximun_iterations);
			}
		}
		//Console.WriteLine("val = " + val + " , st = " + st + " , it = " + i);
		return st;
	}

	/**  Try steffensen acceleration (if error, return -1) 
	 @param f function 
	 @param d first derivative 
	 @param x independent variable 
	 @param val value to add (0 if pure root search) 
	 @param iterations number of iterations 
	 @return  root x value 
	*/
	public final double TrySteffensenAcceleration(function f, derivative d, double x, double val, double iterations)
	{
		double gx = f.invoke(x) - val;
		int i = 0;
		double xi = -1;
		double xiMin1 = -1;
		double xiMin2 = -1;

		while (Math.abs(gx) > epsilon)
		{
			if (d.invoke(x) == 0)
			{
				return 0;
			}
			if (i < 4 || i % 2 == 0 || xi - 2 * xiMin1 + xiMin2 == 0)
			{
				x = x - gx / d.invoke(x);
			}
			else
			{
				x = GetAitkenIteration(xiMin2, xiMin1, xi);
			}
			gx = f.invoke(x) - val;
			xiMin2 = xiMin1;
			xiMin1 = xi;
			xi = x;
			i++;
			if (i > maxIterations)
			{
				throw new RuntimeException(Strings.Exceed_the_maximun_iterations);
			}
		}
		iterations = i;
		return x;
	}

	//endregion

	//region Search Root with derivatives calculation

	/**  Search one root in an interval 
	 @param X independent variable time series 
	 @param Y dependent variable time series 
	 @param min min value of interval 
	 @param max max value of interval 
	 @param val value to search 
	 @return  x value 
	*/
	public final double SearchOneRoot(ArrayList<Double> X, ArrayList<Double> Y, double min, double max, double val)
	{
		if (X.size() != Y.size())
		{
			throw new IllegalArgumentException(Strings.X_must_have_same_size_as_Y);
		}

		Splines sp = new Splines(X, Y, true);
		fCol = new HashMap<String, Double>();
		//cargar la f con los sp.Interpolar

		Derivatives de = new Derivatives();
		//fdxCol = de.GetDerivative(X,Y,1); cargar diccionario
		//fdx2Col = de.GetDerivative(X,Y,2); cargar diccionario
		return SteffensenAccOneRoot((double x) -> f(x), (double x) -> fdx(x), (double x) -> fdx2(x), min, max, val);
	}

	//endregion

	//region Fill Percentiles


	/**  Set all values within an interval 
	 @param values dictionary containing values (return parameter) 
	 @param f function 
	 @param min min value of independent variable 
	 @param max max value of independent variable 
	 @param fmin min value of dependent variable 
	 @param fmax max value of dependent variable 
	*/
	public final int SetAllValues(HashMap<Integer, Double> values, function f, double min, double max, double fmin, double fmax, boolean bisection)
	{
		if (bisection)
		{
			return SetAllValuesBisection(values, (double x) -> f(x), min, max, fmin, fmax);
		}
		else
		{
			return SetAllValuesMontecarlo(values, (double x) -> f(x), min, max, fmin, fmax);
		}
	}

	private int SetAllValuesMontecarlo(HashMap<Integer, Double> values, function f, double min, double max, double fmin, double fmax)
	{
		int its = 0;
		tangible.RefObject<Integer> tempRef_its = new tangible.RefObject<Integer>(its);
		HashMap<Double, Double> rs = RangeSearch((double x) -> f(x), fmin, fmax, min, max, 0.1, epsilon, tempRef_its, maxIterations);
	its = tempRef_its.argValue;
		for (double key : rs.keySet())
		{
			values.put(GetKey(key), rs.get(key));
		}
		return its;
	}

	private int SetAllValuesBisection(HashMap<Integer, Double> values, function f, double min, double max, double fmin, double fmax)
	{
		int it = 0;
		int its = 0;
		if (max - min <= 0.1)
		{
			values.put(GetKey(min), fmin);
			values.put(GetKey(max), fmax);
		}
		else if (fmax - fmin < 1)
		{
			for (int i = GetKey(min);i <= GetKey(max);i++)
			{
				values.put(i, fmax);
			}
		}
		else
		{
			double fmid = (fmin + fmax) / 2.0;
			double mid = Math.round((f.invoke(fmid) * 100) * Math.pow(10, 1)) / Math.pow(10, 1);
			if (mid <= min || mid >= max)
			{
				mid = (min + max) / 2.0;
				tangible.RefObject<Integer> tempRef_it = new tangible.RefObject<Integer>(it);
				fmid = MonotoneBisection((double x) -> f(x), true, fmin, fmax, mid / 100.0, 0.001, tempRef_it, 100);
			it = tempRef_it.argValue;
				its += it;
			}
			its += SetAllValuesBisection(values, (double x) -> f(x), min, mid, fmin, fmid);
			its += SetAllValuesBisection(values, (double x) -> f(x), mid, max, fmid, fmax);
		}
		return its;
	}

	private int GetKey(double val)
	{
		return (int)java.lang.Math.round(Math.round(val * Math.pow(10, 1)) / Math.pow(10, 1) * 10);
	}


	//endregion

	//region Fourier Methods

	/** 
	 <p>Devuelve el punto de fourier con el que el metodo de Newton Raphson converge muy rapidamente</p>
	 <p>Fourier Convergence Conditions value. Returns zero if does not fit conditions</p>
	 <p>1) f(a)*f(b) less than 0   2) f" greater than 0 en [a,b]  val is for non roots (otherwise 0)</para>
	 
	 @param f funcion mediante la cual se calcula el valor de dicha funcion
	 @param d2 funcion mediante la cual se calcula el valor de la derivada segunda de dicha funcion
	 @param a Valor minimo del intervalo en el cual busca la raiz o cero
	 @param b Valor maximo del intervalo en el cual busca la raiz o cero
	 @param val Valor para el cual se desea buscar la inversa de la función <paramref name="f"/>
	 @return  the calculated fourier value 
	*/
	private double GetValueInFourierInterval(function f, derivative2 d2, double a, double b, double val)
	{
		double posValue = 0;
		double negValue = 0;
		double ga = f.invoke(a) - val;
		double gb = f.invoke(b) - val;
		if (ga * gb >= 0)
		{
			return 0;
		}
		if (ga > 0)
		{
			posValue = a;
			negValue = b;
		}
		else
		{
			posValue = b;
			negValue = a;
		}

		//Second derivative constant in the interval (warning: sign is not verified on intermediate values)
		if (d2.invoke(a) > 0 && d2.invoke(b) > 0)
		{
			return posValue;
		}
		else
		{
			return negValue;
		}
	}

	/** 
	 <p>Encuentra el punto de fourier con el que el metodo de Newton Raphson converge muy rapidamente</p>
	 <p>Fourier Convergence Conditions value for functions with one root on interval</p>
	 
	 @param f funcion mediante la cual se calcula el valor de dicha funcion
	 @param d2 funcion mediante la cual se calcula el valor de la derivada segunda de dicha funcion
	 @param min Valor minimo del intervalo en el cual busca la raiz o cero
	 @param max Valor maximo del intervalo en el cual busca la raiz o cero
	 @param val Valor para el cual se desea buscar la inversa de la función <paramref name="f"/>
	 Si no existe una raiz, entonces genera error
	 @return  the calculated fourier value 
	*/
	public final double GetFourierVal(function f, derivative2 d2, double min, double max, double val)
	{
		double a = randGen.NextDouble(min, max);
		double b = randGen.NextDouble(min, max);
		int i = 0;
		if (i > maxIterations - 1)
		{
			a = min;
			b = max;
		}
		double ga = f.invoke(a) - val;
		double gb = f.invoke(b) - val;
		if (ga * gb < 0)
		{
			return GetValueInFourierInterval((double x) -> f(x), d2, a, b, val);
		}

		while (true)
		{
			a = randGen.NextDouble(min, max);
			b = randGen.NextDouble(min, max);
			ga = f.invoke(a) - val;
			gb = f.invoke(b) - val;
			if (ga * gb < 0)
			{
				return GetValueInFourierInterval((double x) -> f(x), d2, a, b, val);
			}
			if (i > maxIterations)
			{
				throw new RuntimeException(Strings.Exceed_the_maximun_iterations);
			}
			i++;
		}
	}

	/** 
	 <p>Encuentra el punto de fourier con el que el metodo de Newton Raphson converge muy rapidamente</p>
	 <p>Fourier Convergence Conditions value for functions with one root on interval</p>
	 
	 @param f funcion mediante la cual se calcula el valor de dicha funcion
	 @param d2 funcion mediante la cual se calcula el valor de la derivada segunda de dicha funcion
	 @param min Valor minimo del intervalo en el cual busca la raiz o cero
	 @param max Valor maximo del intervalo en el cual busca la raiz o cero
	 @param val Valor para el cual se desea buscar la inversa de la función <paramref name="f"/>
	 Si no existe una raiz, entonces genera error
	 @return  the calculated fourier value 
	*/
	public final double GetFourierValue(function f, derivative2 d2, double min, double max, double val)
	{
		return GetFourierValue((double x) -> f(x), d2, min, max, val, 0);
	}

	private double GetFourierValue(function f, derivative2 d2, double min, double max, double val, int it)
	{
		it++;
		if (it > maxIterations)
		{
			//Console.WriteLine("Obtener Fourier Value: Mas de 30 iteraciones");
			return GetValueInFourierInterval((double x) -> f(x), d2, min, max, val);
		}
		double gMin = f.invoke(min) - val;
		double gMax = f.invoke(max) - val;
		double inter = min + (max - min) / 2.0;
		double gInter = f.invoke(inter) - val;
		if (Math.abs(gInter) < 0.02)
		{
			return inter;
		}
		if (gMin * gInter < 0)
		{
			return GetFourierValue((double x) -> f(x), d2, min, inter, val, it);
		}
		else
		{
			return GetFourierValue((double x) -> f(x), d2, inter, max, val, it);
		}
	}

	//endregion

	//region Auxiliar Methods

	/**  Quadratic function, just for performance 
	 @param val value to square 
	 @return  squared value 
	*/
	public final double quad(double val)
	{
		return val * val;
	}

	//endregion

	//region Class Cell

	public static class Range
	{
		public double xIni;
		public double xEnd;
		public double yIni;
		public double yEnd;
		public int nValues;

		public Range(double xIni, double xEnd, double yIni, double yEnd, int nValues)
		{
			this.xIni = xIni;
			this.xEnd = xEnd;
			this.yIni = yIni;
			this.yEnd = yEnd;
			this.nValues = nValues;
		}
	}

	//endregion

}