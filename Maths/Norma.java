package NumCalc;

import ClassicalStat.*;
import java.util.*;

//region Imports


//endregion



/**  Data normalization class 
*/
public class Norma
{

	//region Fields

	private double[] datos;
	private double[] norm;
	private double min;
	private double max;
	private double normMax;
	private double rango;
	private int maxClasses;
	private NormType type = NormType.values()[0];
	private Statistics stat;
	private List<Double> serie;

	//endregion

	//region Constructors

	 /**  Constructor por defecto 
	 */
	public Norma()
	{
		type = NormType.none;
		normMax = 1.0;
		stat = new Statistics();
	}

	/**  Constructor con parametro 
	 @param datos data to normalize 
	*/
	public Norma(double[] datos)
	{
		this.datos = datos;
		type = NormType.none;
		normMax = 1.0;
	}

	/**  Constructor con parametro 
	 @param datos data to normalize 
	*/
	public Norma(ArrayList<Double> datos)
	{
		this.datos = tangible.DoubleLists.toArray(datos);
		type = NormType.none;
		normMax = 1.0;
	}

	/**  Constructor con parametro 
	 @param datos data to normalize 
	*/
	public Norma(List datos)
	{
		this.datos = new double[datos.size()];
		for (int i = 0;i < datos.size();i++)
		{
			this.datos[i] = (double)datos.get(i);
		}
		normMax = 1.0;
	}

	public Norma(List<Double> datos)
	{
		this.serie = datos;
		this.datos = new double[datos.size()];
		for (int i = 0; i < datos.size(); i++)
		{
			this.datos[i] = (double)datos.get(i);
		}
		normMax = 1.0;
	}

	//endregion

	//region Properties

	/**  datos crudos 
	*/
	public final double[] getDatos()
	{
		return datos;
	}
	public final void setDatos(double[] value)
	{
		datos = value;
	}

	/**  datos normalizados 
	*/
	public final double[] getNorm()
	{
		return norm;
	}

	/**  Tipo de normalizacion 
	*/
	public final NormType getType()
	{
		return type;
	}

	/**  Minimo 
	*/
	public final double getMin()
	{
		return min;
	}

	/**  Maximo 
	*/
	public final double getMax()
	{
		return max;
	}

	/**  Nuevo Maximo 
	*/
	public final double getNormMax()
	{
		return normMax;
	}
	public final void setNormMax(double value)
	{
		normMax = value;
	}

	/**  Rango 
	*/
	public final double getRango()
	{
		return rango;
	}

	/**  MaxClasses 
	*/
	public final int getMaxClasses()
	{
		return maxClasses;
	}
	public final void setMaxClasses(int value)
	{
		maxClasses = value;
	}

	//endregion

	//region Public Methods

	/**  Normalize data 
	 @param type type of normalization 
	*/
	public final void Normalize(NormType type)
	{
		this.type = type;
		switch (type)
		{
			case none:
				this.norm = this.datos;
				break;
			case minMax:
				NormalizeMinMax();
				break;
			case log:
				NormalizeLog();
				break;
			case logMinMax:
				NormalizeLogMinMax();
				break;
			case min:
				NormalizeMin();
				break;
			case max:
				NormalizeMax();
				break;
			case hash:
				NormalizeHash();
				break;
			case hashMinMax:
				NormalizeHashMinMax();
				break;
			case scale:
				NormalizeScale();
				break;
			case signMax:
				NormalizeSignMax();
				break;

		}
	}

	/**  Get Data 
	 @param index index of datum 
	 @param normalized if normalized 
	 @return  value of requested position 
	*/
	public final double GetValue(int index, boolean normalized)
	{
		if (!normalized)
		{
			return datos[index];
		}
		else
		{
			return norm[index];
		}
	}

	//endregion

	//region Private Methods

	private void SetMinMax()
	{
		min = Double.MAX_VALUE;
		max = -Double.MAX_VALUE;
		for (int i = 0;i < datos.length;i++)
		{
			if (datos[i] < min)
			{
				min = datos[i];
			}
			if (datos[i] > max)
			{
				max = datos[i];
			}
		}
		rango = max - min;
	}

	private void NormalizeMinMax()
	{
		SetMinMax();
		if (rango == 0)
		{
			norm = datos;
			return;
		}
		norm = new double[datos.length];
		for (int i = 0;i < datos.length;i++)
		{
			norm[i] = ((datos[i] - min) / rango) * normMax;
		}
	}

	private void NormalizeMin()
	{
		SetMinMax();
		norm = new double[datos.length];
		for (int i = 0; i < datos.length; i++)
		{
			norm[i] = (datos[i] - min);
		}
	}

	private void NormalizeMax()
	{
		SetMinMax();
		norm = new double[datos.length];
		for (int i = 0;i < datos.length;i++)
		{
			norm[i] = datos[i] / max;
		}
	}

	private void NormalizeLog()
	{
		norm = new double[datos.length];
		for (int i = 0;i < datos.length;i++)
		{
		   norm[i] = Math.log(datos[i]);
		}
	}

	private void NormalizeLogMinMax()
	{
		SetMinMax();
		norm = new double[datos.length];
		for (int i = 0;i < datos.length;i++)
		{
			norm[i] = Math.log((datos[i] - min) / rango);
		}
	}

	private void NormalizeHash()
	{
		double mean = stat.Mean(datos);
		double stDev = stat.StDev(datos);
		norm = new double[datos.length];
		for (int i = 0;i < datos.length;i++)
		{
			norm[i] = (datos[i] - mean) / Math.pow(stDev,2);
		}
	}

	private void NormalizeHashMinMax()
	{
		NormalizeHash();
		double min = Double.MAX_VALUE;
		double max = -Double.MAX_VALUE;
		for (int i = 0;i < norm.length;i++)
		{
			if (norm[i] < min)
			{
				min = norm[i];
			}
			if (norm[i] > max)
			{
				max = norm[i];
			}
		}
		rango = max - min;
		for (int i = 0;i < norm.length;i++)
		{
			norm[i] = (norm[i] - min) / rango;
		}
	}

	private void NormalizeScale()
	{

		SetMinMax();
		if (rango < maxClasses)
		{
			norm = datos;
			return;
		}
		double[] normMinMax = new double[datos.length];
		for (int i = 0;i < datos.length;i++)
		{
			normMinMax[i] = Math.ceil(((datos[i] - min) / rango) * 100);
		}

		norm = new double[normMinMax.length];
		int width = (int)(Math.ceil((double)datos.length / (double)maxClasses));
		int iniClass = 0;
		int nClass = 0;
		double totClass = 0.0;
		double mean;
		for (int i = 0;i < normMinMax.length;i++)
		{
			totClass += normMinMax[i];
			nClass++;
			if (nClass == width || i == normMinMax.length - 1)
			{
				mean = totClass / nClass;
				for (int j = iniClass;j <= i;j++)
				{
					norm[j] = mean;
				}
				nClass = 0;
				totClass = 0.0;
				iniClass = i + 1;
			}
		}
	}

	private void NormalizeSignMax()
	{
		SetMinMax();
		if (Math.abs(min) > Math.abs(max))
		{
			max = Math.abs(min);
		}
		min = 0;
		rango = max;
		norm = new double[datos.length];
		for (int i = 0; i < datos.length; i++)
		{
			norm[i] = datos[i] / max;
		}


	}

	public final double Normalize(double value)
	{
	   return ((value - min) / rango) * normMax;
	}

	public final double UnNormalize(double value)
	{
		return (value / normMax) * (max - min) + min;
	}

	/**  Get scaled value 
	 @param unSc the scale to apply 
	 @return  scaled value 
	*/
	public final double GetScaled(double unSc)
	{
		return Math.ceil((unSc - min) / rango) * 100;
	}

	/**  Get unscaled value 
	 @param sc the scale to apply 
	 @return  unscaled value 
	*/
	public final double GetUnScaled(double sc)
	{
		return ((sc / 100.0) * rango) + min;
	}

	//endregion

	//region Enums

	/**  Type of normalization 
	*/
	public enum NormType
	{
		/**  No normalization 
		*/
		none,
		/**  Min-Max normalization (from 0 to 1) 
		*/
		minMax,
		/**  Logaritmized normalization 
		*/
		log,
		/**  Logaritmized and then Min-Max 
		*/
		logMinMax,
		/**  Only min normalization 
		*/
		min,
		/**  Only max normalization 
		*/
		max,
		/**  Hash normalization, optimizing hash level 
		*/
		hash,
		/**  Hash normalization and then  min-max normalization 
		*/
		hashMinMax,
		/**  Scale normalization 
		*/
		scale,
		/**  Max normalization with sign 
		*/
		signMax;

		public int getValue()
		{
			return this.ordinal();
		}

		public static NormType forValue(int value)
		{
			return values()[value];
		}
	}

	//endregion
}