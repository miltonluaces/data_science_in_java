package DeepLearning;

import RedesNeuronales.Properties.*;

//region Imports


//endregion

/**  Class that represent a process element unit (Neuron) 
*/
public class Neurona
{

	//region Fields

	private int red; // identificador  red
	private int cant; // cantidad dendritas
	private int capa; // capa en la red
	private int indice; // indice en la capa
	private double n; // coeficiente de entrenamiento
	private double m; // factor momento
	private double[][] matriz;
	private NeuronaCollection conexiones;
	private int cantConexiones;
	private Funcion funcionActivacion = Funcion.sigmoidal;

	//endregion

	//region Constructor

	/**  Constructor 
	 @param red number of net 
	 @param capa number of layer 
	 @param indice number of index 
	 @param cant quantity of neurons 
	*/
	public Neurona(int red, int capa, int indice, int cant)
	{
	   this.red = red;
	   this.capa = capa;
	   this.indice = indice;
	   this.cant = cant;
	   this.n = 1.0;
	   this.m = 1.0;
	   conexiones = new NeuronaCollection();
	   cantConexiones = cant;

	   matriz = new double[cant][4];
	   for (int i = 0;i < cant;i++)
	   {
		   for (int j = 0;j < 3;j++)
		   {
				 matriz[i][j] = 0.00;
		   }
		   matriz[i][3] = red * 100000 + capa * 10000 + indice * 100 + i; //id basado en red.capa.indice.dendrita
	   }

	   if (capa != 0)
	   {
		   matriz[0][0] = 1.00;
	   } //neurona marcadora de tendencia
	}

	//endregion

	//region Properties

	/**  Number of the net 
	*/
	public final int getRed()
	{
		return red;
	}
	public final void setRed(int value)
	{
		red = value;
	}

	/**  Number of the layer 
	*/
	public final int getCapa()
	{
		return capa;
	}
	public final void setCapa(int value)
	{
		capa = value;
	}

	/**  Index of the neuron 
	*/
	public final int getIndice()
	{
		return indice;
	}
	public final void setIndice(int value)
	{
		indice = value;
	}

	/**  trainning coefficent 
	*/
	public final double getM()
	{
		return m;
	}
	public final void setM(double value)
	{
		m = value;
	}

	/**  moment factor 
	*/
	public final double getN()
	{
		return n;
	}
	public final void setN(double value)
	{
		n = value;
	}

	/**  number of connections 
	*/
	public final int getCantConexiones()
	{
		return cantConexiones;
	}

	public final Funcion getFuncionActivacion()
	{
		return funcionActivacion;
	}
	public final void setFuncionActivacion(Funcion value)
	{
		funcionActivacion = value;
	}

	//endregion

	//region Setters and Getters

	/**  Set value 
	 @param unX the value 
	 @param unIndice index in the array 
	*/
	public final void setX(double unX, int unIndice)
	{
		matriz[unIndice][0] = unX;
	}

	/**  Get value 
	 @param unIndice index in the array 
	  the value 
	*/
	public final double getX(int unIndice)
	{
		return matriz[unIndice][0];
	}


	/**  Set weight 
	 @param unW the weight 
	 @param unIndice index in the array 
	*/
	public final void setW(double unW, int unIndice)
	{
		matriz[unIndice][1] = unW;
	}

	/**  Get weight 
	 @param unIndice index in the array 
	  the weight 
	*/
	public final double getW(int unIndice)
	{
		return matriz[unIndice][1];
	}

	/**  Set previous weight 
	 @param unWAnt the weight 
	 @param unIndice index in the array 
	*/
	public final void setWAnt(double unWAnt, int unIndice)
	{
		matriz[unIndice][2] = unWAnt;
	}

	/**  Get previous weight 
	 @param unIndice index in the array 
	  the weight 
	*/
	public final double getWAnt(int unIndice)
	{
		return matriz[unIndice][2];
	}

	/**  Set ID 
	 @param unID the ID 
	 @param unIndice index in the array 
	*/
	public final void setID(long unID, int unIndice)
	{
		matriz[unIndice][3] = unID;
	}

	/**  Get ID 
	 @param unIndice index in the array 
	 @return  the ID 
	*/
	private long getID(int unIndice)
	{
		return (long)matriz[unIndice][3];
	}

	/**  Get Weight Matrix  
	 @return  the matrix 
	*/
	public final double[][] getMatriz()
	{
		return matriz;
	}

	/**  Add a new forward connection 
	 @param ne the target Neuron 
	*/
	public final void agregarConexion(Neurona ne)
	{
		conexiones.agregar(ne);
		cantConexiones++;
	}

	/**  Get all connections 
	 @return  a list of collections 
	*/
	public final NeuronaCollection getConexiones()
	{
		 return conexiones;
	}

	//endregion        

	//region Entrenamiento y Proceso Internos

	public final double axon()
	{
		 return GetFuncionActivacion(net());
	}

	public final void transmitir()
	{
		double y;
		if (capa == 0)
		{
			y = matriz[0][0];
		} //unica entrada
		else
		{
			y = axon();
		}
		for (Neurona ne : conexiones)
		{
			ne.setX(y, indice);
		}
	}

	public final void ajustarPesos(double d)
	{
		double pesoActual;
		for (int i = 0;i < cant;i++)
		{
			pesoActual = getW(i);
			setW(getW(i) + n * getDelta(d) * getX(i) + m * (getW(i) - getWAnt(i)), i); //n�x + momento
			setWAnt(pesoActual, i); //ahora guardo el peso anterior
		}
	}

	public final void ajustarPesos(double[] ds)
	{
		double pesoActual;
		for (int i = 0;i < cant;i++)
		{
			pesoActual = getW(i);
			setW(getW(i) + n * getDelta(ds) * getX(i) + m * (getW(i) - getWAnt(i)), i); //n�x + momento
			setWAnt(pesoActual, i); //ahora guardo el peso anterior
		}
	}

	//endregion

	//region Private Methods

	private double net()
	{
		double net = 0;
		for (int i = 0;i < cant;i++)
		{
			net = net + matriz[i][0] * matriz[i][1];
		}
		return net;
	}

	private double GetFuncionActivacion(double net)
	{
		switch (funcionActivacion)
		{
			case sigmoidal:
				double MenosNet = - net;
				double exp = Math.exp(MenosNet);
				double res = 1 / (1 + exp);
				return res;
			case tangenteHiperbolica:
				return Math.tanh(net);
			default:
				throw new UnsupportedOperationException(Strings.Function_not_implemented);
		}
	}

	private double getDelta(double d)
	{
		double delta = 0.00;
		double y = axon();

		if (capa == 2)
		{
			delta = y * (1 - y) * (d - y); // y'(d-y)
		}

		else if (capa == 1)
		{
			for (Neurona ne : conexiones)
			{
				delta = y * (1 - y) * (ne.getW(indice) * ne.getDelta(d));
			}
		}
		//Console.WriteLine("delta = " + delta);
		return delta;
	}

	private double getDelta(double[] ds)
	{
		double delta = 0.00;
		double y = axon();

		if (capa == 2)
		{
			delta = y * (1 - y) * (ds[indice] - y); // y'(d-y)
		}

		else if (capa == 1)
		{
			for (Neurona neu : this.conexiones)
			{
				delta = delta + y * (1 - y) * (neu.getW(indice) * neu.getDelta(ds));
			}
		}
		return delta;
	}

	//endregion

	//region Public Enums

	public enum Funcion
	{
		sigmoidal,
		tangenteHiperbolica;

		public int getValue()
		{
			return this.ordinal();
		}

		public static Funcion forValue(int value)
		{
			return values()[value];
		}
	}

	//endregion

}