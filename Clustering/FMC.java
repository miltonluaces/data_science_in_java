package Clustering;

import NumCalc.*;

/**  Algorithmic repository for Fuzzy Mean Clustering  
*/
public class FMC
{

	//region Fields

	private Functions func;
	private double[][] VectorDatos;
	private DatoFuzzy[] DatoFuzzy;
	private int cantDatos;
	private int cantClusters;
	private int cantElementos;
	private Cluster[] Cluster;
	private double difCentros;
	private String path;
	private double q; //fuzziness

	//endregion

	//region Constructor

	/**  Constructor 
	 @param unD data matrix 
	 @param unC quantity of data 
	 @param unQ q 
	*/
	public FMC(double[][] unD, int unC, int unQ)
	{
		func = new Functions();
		VectorDatos = unD;
		cantDatos = VectorDatos.length;
		DatoFuzzy = new DatoFuzzy[cantDatos];
		for (int i = 0;i < cantDatos;i++)
		{
			DatoFuzzy[i] = new DatoFuzzy(VectorDatos[i], unC);
		}
		cantElementos = VectorDatos[0].length;
		cantClusters = unC;
		q = unQ;
		difCentros = Double.MAX_VALUE;
		path = "C:\\fuzzyMeans.txt";
		DatoFuzzy centro;
		Cluster = new Cluster[cantClusters];
		for (int i = 0;i < cantClusters;i++)
		{
			centro = new DatoFuzzy(cantElementos);
			Cluster[i] = new Cluster(centro);
		}
	}

	//endregion

	//region Setters and Getters

	/**  Set path for files 
	 @param unP the path to set 
	*/
	public final void setPath(String unP)
	{
		path = unP;
	}

	/**  Get difference between centres 
	 @return  differece between centres 
	*/
	public final double getDifCentros()
	{
		return difCentros;
	}

	//endregion

	//region Metodos de Calculo

	/**  Fuzzy means first step 
	*/
	public final void calcularCentros()
	{
		//ajuste global para todos los centros
		double num = 0.0;
		double den = 0.0; //         \� K      q          [df,K = cant datos]
		DatoFuzzy centro; //         /_ df=1  M ik  U ke   [e,E = cant elementos]
		double mq; // C ie = �������������������   [U ke = elemento e del dato df]
		for (int i = 0;i < cantClusters;i++)
		{ //            \� K      q       [i,I = cant clusters]
											   //            /_ df=1  M ik      [M ik = membership del dato df al cluster i]
			centro = new DatoFuzzy(cantElementos);
			for (int e = 0;e < cantElementos;e++)
			{
				for (int k = 0;k < cantDatos;k++)
				{
					mq = Math.pow(DatoFuzzy[k].GetM(i), q);
					num += mq * DatoFuzzy[k][e];
					den += mq;
				}
				centro[e] = num / den;
			}
			Cluster[i].Centro = centro;
		}
	}

	/**  Fuzzy means second step 
	*/
	public final void asignarMembership()
	{
		double distACentro_i; //                 1
		double sumatoria = 0f; // M ik =  ��������������������
		double mik; //           /_ j=1  [D ik/ D jk]
		//double totMik = 0.0;
		for (int i = 0;i < cantClusters;i++)
		{
			for (int k = 0;k < cantDatos;k++)
			{
				distACentro_i = func.DistEuclid(DatoFuzzy[k].Valores,Cluster[i].Centro.Valores);
				sumatoria = 0.0;
				for (int j = 0;j < cantClusters;j++)
				{
					{
					/*if(j != i)*/   sumatoria += Math.pow((distACentro_i / func.DistEuclid(DatoFuzzy[k].Valores,Cluster[j].Centro.Valores)),2 / (q - 1));
					}
				}
				mik = 1.00 / sumatoria;
				//if(df==0) totMik += mik;
				//System.out.println(""+mik);
				//if(df==0) System.out.println("Tot "+df+" : "+ totMik);
				DatoFuzzy[k].SetM(mik,i);
			}
		}
	}

	/**  Main calculation method 
	*/
	public final void clustering()
	{
		while (difCentros > 0)
		{
			calcularCentros();
			asignarMembership();
		}
	}

	//endregion

	//region Log Generation 

	/**  Generate log 
	*/
	public final void clusteringLog()
	{
		try
		{
//C# TO JAVA CONVERTER WARNING: The java.io.OutputStreamWriter constructor does not accept all the arguments passed to the System.IO.StreamWriter constructor:
//ORIGINAL LINE: StreamWriter so = new StreamWriter(path,false);
			java.io.OutputStreamWriter so = new java.io.OutputStreamWriter(path);
			so.write("CLUSTERING - Fuzzy-MEANS : " + cantDatos + " Datos con " + cantElementos + " elementos en " + cantClusters + " clusters." + System.lineSeparator());
			so.write("------------------------------------------------------------------------------------------------------" + System.lineSeparator());
			//while(difCentros > 0.5)
			for (int j = 0;j < 4;j++)
			{
				so.write("" + System.lineSeparator());
				so.write("ITERACION " + j + ":" + System.lineSeparator());
				so.write("" + System.lineSeparator());
				calcularCentros();
				asignarMembership();
				for (int k = 0;k < cantDatos;k++)
				{
					so.write(String.valueOf(DatoFuzzy[k]) + System.lineSeparator());
				}
			}
			so.close();
			System.Diagnostics.Debug.WriteLine("Completado archivo " + path);
		}
		catch (java.lang.Exception e)
		{
			System.Diagnostics.Debug.WriteLine("Error");
		}
	}

	//endregion

	//region ToString Override

	/**  ToString override  
	 @return  the formatted string 
	*/
	@Override
	public String toString()
	{
		String ret = String.format(Properties.Strings.Movimientos_centros, difCentros);
		for (int i = 0; i < cantClusters; i++)
		{
			ret = ret + "\n" + Properties.Strings.CLUSTER + (i + 1) + "\n" + Cluster[i] + "\n";
		}
		return ret;
	}

	//endregion
}