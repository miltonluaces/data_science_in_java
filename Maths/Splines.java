package NumCalc;

//region Imports

import java.util.*;

//endregion

/**  Linear and cubic b-splines Class 
*/
public class Splines
{

	//region Fields

	private ArrayList<Point> points;
	private double[] y2;
	private ArrayList<Point> rectas;
	private boolean linear;

	//endregion

	//region Constructors

	/**  Default constructor 
	*/
	public Splines()
	{
	}

	/** 
	 Constructor que recibe un arraylist de elementos Point (x,y)
	 
	 @param points list of points 
	 @param linear if it is linear spline 
	*/

	public Splines(ArrayList<Point> points, boolean linear)
	{
		this.points = points;
		this.linear = linear;
		if (linear)
		{
			PrecalculoLinear();
		}
		else
		{
			PrecalculoCubica();
		}
	}

	/** 
	 Constructor que recibe vectores de x e y
	 
	 @param x independent variable 
	 @param y dependent variable 
	 @param linear if it is linear spline 
	*/
	public Splines(ArrayList<Double> x, ArrayList<Double> y, boolean linear)
	{
		if (x.size() != y.size())
		{
			throw new RuntimeException(Strings.los_arrays_deben_tener_el_mismo_tama�o);
		}
		this.points = new ArrayList<Point>();
		for (int i = 0;i < x.size();i++)
		{
			points.add(new Point((double)x.get(i), (double)y.get(i)));
		}
		this.linear = linear;
		if (linear)
		{
			PrecalculoLinear();
		}
		else
		{
			PrecalculoCubica();
		}
	}

	//endregion

	//region Public Methods

	/** 
	 Metodo Principal de interpolacion
	 Para cada valor de x, rectas[x] tiene los coeficientes a y b de la recta (en x e y)
	 
	 @param x independent variable for interpolation 
	 @return  interpolation 
	*/
	public final double Interpolar(double x)
	{
		double interp = 0.0;
		if (linear)
		{
			if (x >= ((Point)points.get(points.size() - 1)).x)
			{
				return ((Point)points.get(points.size() - 1)).y;
			}
			else
			{
				return ((Point)rectas.get((int)x)).x * x + ((Point)rectas.get((int)x)).y;
			}
		}
		else
		{
			int n = this.points.size();
			int klo = 0;
			int khi = n - 1;

			//biseccion
			while (khi - klo > 1)
			{
//C# TO JAVA CONVERTER WARNING: The right shift operator was not replaced by Java's logical right shift operator since the left operand was not confirmed to be of an unsigned type, but you should review whether the logical right shift operator (>>>) is more appropriate:
				int k = (khi + klo) >> 1;
				if (GetX(k) > x)
				{
					khi = k;
				}
				else
				{
					klo = k;
				}
			}

			double h = GetX(khi) - GetX(klo);
			double a = (GetX(khi) - x) / h;
			double b = (x - GetX(klo)) / h;

			//Evaluacion del spline cubico
			interp = a * GetY(klo) + b * GetY(khi) + ((a * a * a - a) * y2[klo] + (b * b * b - b) * y2[khi]) * (h * h) / 6.0;
		}
		return interp;
	}

	/** 
	 Interpolacion puntual de un punto en un intervalo
	 
	 @param x independent variable for interpolation 
	 @param x1 first defining point (x axis value) 
	 @param y1 first defining point (y axis value) 
	 @param x2 second defining point (x axis value) 
	 @param y2 second defining point (y axis value)
	 @return  functional value for x 
	*/
	public final double Interpolar(double x, double x1, double y1, double x2, double y2)
	{
		double a = (y2 - y1) / (x2 - x1);
		double b = y1 - (x1 * (y2 - y1) / (x2 - x1));
		return a * x + b;
	}

	/**  Obtener las n inversas de los splines 
	 x = (y - b) / a
	 @param y dependent variable 
	 @return  indepenedent variables 
	*/
	public final ArrayList<Double> Inversas(double y)
	{
		ArrayList<Double> inversas = new ArrayList<Double>();
		Point recta = new Point();
		double x;
		for (int i = 0;i < rectas.size();i++)
		{
			recta = (Point)rectas.get(i);
			if (recta.x == 0)
			{ // si tiene derivada nula
				if (y == recta.y)
				{
				}
			}
			else
			{
				x = (y - recta.y) / recta.x; //si intercepta el liston
				if (x >= i && x < i + 1)
				{
					inversas.add(x);
				}
			}
		}
		return inversas;
	}

	/**  Positive intervals calculation 
	 @param y dependent variable 
	 @return  list of points that define the intervals 
	*/
	public final ArrayList<Point> IntervalosPositivos(double y)
	{
		ArrayList<Point> intervalos = new ArrayList<Point>();
		Point recta = new Point();
		double x;
		double delta = 0.00001;
		boolean meseta = false;
		Point intervalo = new Point();
		for (int i = 0;i < rectas.size();i++)
		{
			recta = (Point)rectas.get(i);
			if (recta.x == 0)
			{ //derivada cero, es meseta
				if (y == recta.y)
				{
					intervalo = new Point(i, i + 1);
					intervalos.add(intervalo.clone());
					meseta = true;
				}
			}
			else
			{
				x = (y - recta.y) / recta.x; //es punto de interseccion
				if (x >= i - delta && x <= (i + 1) + delta)
				{ //es valido
					if (recta.x > 0)
					{ //derivada positiva
						if (meseta)
						{
							meseta = false;
						}
						intervalo = new Point();
						intervalo.x = x; //el punto de interseccion es apertura de intervalo
					}
					else
					{ //derivada negativa
						if (meseta)
						{
							meseta = false;
						}
						else
						{
							intervalo.y = x; //el punto de interseccion es cierre de intervalo
							intervalos.add(intervalo.clone());
						}
					}
				}
			}
		}
		if (intervalos.isEmpty())
		{
			return intervalos;
		}

		//join de intervalos
		ArrayList<Point> intervalosJoin = new ArrayList<Point>();
		Point act = new Point();
		Point next = new Point();
		act = (Point)intervalos.get(0);
		double a = act.x;
		double b;
		meseta = false;
		for (int i = 0;i < intervalos.size();i++)
		{
			act = (Point)intervalos.get(i);
			if (i < intervalos.size() - 1)
			{
				next = (Point)intervalos.get(i + 1);
				if (Math.abs(act.y - next.x) < delta)
				{
					if (!meseta)
					{
						a = act.x;
						meseta = true;
					}
					continue;
				}
			}
			if (!meseta)
			{
				a = act.x;
			}
			b = act.y;
			intervalosJoin.add(new Point(a, b));
			meseta = false;
		}
		return intervalosJoin;
	}

	/**  Agregar un valor 
	 @param x independent variable 
	 @param y dependent variable 
	*/
	public final void Add(double x, double y)
	{
		Point p = new Point(x, y);
		points.add(p.clone());
	}

	/** 
	 Reset Vector de valores
	*/
	public final void Clear()
	{
		this.points.clear();
	}

	//endregion

	//region Private Methods

	private void PrecalculoLinear()
	{
		rectas = new ArrayList<Point>();
		Point recta = new Point();
		double a;
		double b;
		for (int i = 0;i < points.size() - 1;i++)
		{
			a = (GetY(i + 1) - GetY(i)) / (GetX(i + 1) - GetX(i));
			b = GetY(i) - (GetX(i) * (GetY(i + 1) - GetY(i)) / (GetX(i + 1) - GetX(i)));
			recta = new Point(a, b);
			rectas.add(recta.clone());
		}
	}

	private void PrecalculoCubica()
	{
		int n = this.points.size();
		double[] u = new double[n];
		this.y2 = new double[n];
		u[0] = 0;
		this.y2[0] = 0;

		//Descomposicion tridiagonal. y2 y u variables temporales.
		for (int i = 1;i < n - 1;++i)
		{
			double wx = GetX(i + 1) - GetX(i - 1);
			double sig = (GetX(i) - GetX(i - 1)) / wx;
			double p = sig * y2[i - 1] + 2.0;

			this.y2[i] = (sig - 1.0) / p;

			double ddydx = (GetY(i + 1) - GetY(i)) / (GetX(i + 1) - GetX(i)) - ((GetY(i) - GetY(i - 1)) / (GetX(i) - GetX(i - 1)));
			u[i] = (6.0 * ddydx / wx - sig * u[i - 1]) / p;
		}
		this.y2[n - 1] = 0;

		// This is the backsubstitution loop of the tridiagonal algorithm
		for (int i = n - 2;i >= 0;--i)
		{
			this.y2[i] = this.y2[i] * this.y2[i + 1] + u[i];
		}
	}

	private double GetX(int index)
	{
		return (double)((Point)points.get(index)).x;
	}

	private double GetY(int index)
	{
		return (double)((Point)points.get(index)).y;
	}

	//endregion

	//region Struct Point

	/**  Point struct (x,y) 
	*/
//C# TO JAVA CONVERTER WARNING: Java does not allow user-defined value types. The behavior of this class will differ from the original:
//ORIGINAL LINE: public struct Point
	public final static class Point
	{

		/**  x axis value 
		*/
		public double x;
		/**  y axis value 
		*/
		public double y;

		/**  Constructor 
		 @param x x axis value 
		 @param y y axis value 
		*/
		public Point()
		{
		}

		public Point(double x, double y)
		{
			this.x = x;
			this.y = y;
		}

		public Point clone()
		{
			Point varCopy = new Point();

			varCopy.x = this.x;
			varCopy.y = this.y;

			return varCopy;
		}
	}

	//endregion
}