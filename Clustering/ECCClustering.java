package Clustering;

import Clustering.Properties.*;
import java.util.*;

//region Imports


//endregion

/**  Descripci�n breve de ECCClustering. 
*/
public class ECCClustering
{

	//region Fields

	private Dimension[] dim;
	private Dato[] datos;
	private KMC kmc;
	private ArrayList resultClusters = null;

	private ArrayList<Dimension> scalDim;
	private ArrayList<Dato> scalDatos;
	private ArrayList<CatDimension> catDim;
	private ArrayList<ArrayList<boolean[]>> catDatos;

	//endregion

	//region Constructor

	/**  Constructor 
	*/
	public ECCClustering()
	{
		scalDim = new ArrayList<Dimension>();
		scalDatos = new ArrayList<Dato>();
		catDim = new ArrayList<CatDimension>();
		catDatos = new ArrayList<ArrayList<boolean[]>>();
	}

	//endregion

	//region Load Data

	/**  Load data for clustering from arrays 
	 @param dimensions array of dimensions 
	 @param datos array of data 
	*/
	public final void LoadData(Dimension[] dimensions, Dato[] datos)
	{
		this.dim = dimensions;
		this.datos = datos;
		Normalize();
	}

	/** 
	 Load data including categoric dimensions
	 
	 @param scalDim scalar dimensions
	 @param scalDatos scalar data
	 @param catDim categoric dimensions
	 @param catDatos categoric data
	*/
	public final void LoadData(ArrayList<Dimension> scalDim, ArrayList<Dato> scalDatos, ArrayList<CatDimension> catDim, ArrayList<ArrayList<boolean[]>> catDatos)
	{
		this.scalDim = scalDim;
		this.scalDatos = scalDatos;
		this.catDim = catDim;
		this.catDatos = catDatos;
	}

	/**  Load data for clustering from file 
	 @param nombreArchivo name (path) of the file 
	*/
	public final void LoadData(String nombreArchivo)
	{
		String linea = "";
		String[] tokens;
		int nDatos, nDim;
		try
		{
			java.io.InputStreamReader arch = new java.io.InputStreamReader(nombreArchivo);
			char separador = '\t';

			//Numero de datos y Numero de dimensiones
			linea = arch.ReadLine();
			tokens = linea.split(java.util.regex.Pattern.quote(separador.toString()), -1);
			nDatos = Integer.parseInt(tokens[0]);
			nDim = Integer.parseInt(tokens[1]);
			dim = new Dimension[nDim];
			datos = new Dato[nDatos];

			//Dimensiones
			Dimension dimension;
			for (int i = 0;i < nDim;i++)
			{
				linea = arch.ReadLine();
				tokens = linea.split(java.util.regex.Pattern.quote(separador.toString()), -1);
				dimension = new Dimension();
				dimension.index = i;
				dimension.nDivisions = Integer.parseInt(tokens[0]);
				dimension.logarithm = (tokens[1].equals("Y"))? true : false;
				dimension.normalize = (tokens[2].equals("Y"))? true : false;
				dimension.weight = Double.parseDouble(tokens[3]);
				dimension.sort = (tokens[4].equals("Y"))? true : false;
				dimension.min = Double.MAX_VALUE;
				dimension.max = -Double.MAX_VALUE;
				dim[i] = dimension;
			}

			//Datos
			String codigo;
			double[] valores;

			for (int i = 0;i < nDatos;i++)
			{
				linea = arch.ReadLine();
				if (linea == null)
				{
					throw new RuntimeException(String.format(Strings.Ultima_linea_0_No_hay_suficientes_lineas_1, i,nDatos));
				}
				tokens = linea.split(java.util.regex.Pattern.quote(separador.toString()), -1);
				codigo = tokens[0];
				valores = new double[nDim];
				for (int j = 0; j < nDim; j++)
				{
					try
					{
						valores[j] = Double.parseDouble(tokens[j + 1]);
					}
					catch (RuntimeException e)
					{
						throw new RuntimeException(String.format(Strings.La_cadena_0_no_tiene_el_formato_correcto_en_el_token_1_, linea, tokens[j + 1]), e);
					}
					if (valores[j] < dim[j].min)
					{
						dim[j].min = valores[j];
					}
					if (valores[j] > dim[j].max)
					{
						dim[j].max = valores[j];
					}
				}
				datos[i] = new Dato(codigo, valores);
			}

			Normalize();

			//Console.WriteLine("Dimensiones:"); for(int i=0;i<nDim;i++) { Console.WriteLine(i+") " +  " - No.Div " + GetDimension(i).nDivisions + " - Log? " + GetDimension(i).logarithm + " - Norm? " +   GetDimension(i).normalize +  " - Peso: " + GetDimension(i).weight + " - Sort? " + GetDimension(i).sort); } 
			//Console.WriteLine("Datos:"); for(int fil=0;fil<nDatos;fil++) { for(int col=0;col<nDim;col++) { Console.Write(valores[fil][col]+ "  ");  }  Console.WriteLine("");  }
		}
		catch (FileNotFoundException fe)
		{
			throw new FileNotFoundException(Strings.No_se_encontro_el_archivo + nombreArchivo, fe);
		}
	}

	private void Normalize()
	{
		for (int i = 0;i < dim.length;i++)
		{
			//no se logaritma ni se normaliza
			if (!dim[i].logarithm && !dim[i].normalize)
			{
				for (Dato dato : datos)
				{
					dato.Norm[i] = dato[i];
				}
			}
				//se logaritma pero no se normaliza
			else if (dim[i].logarithm && !dim[i].normalize)
			{
				for (Dato dato : datos)
				{
					dato.Norm[i] = Math.log(dato[i]);
				}
			}
				//no se logaritma pero se normaliza
			else if (!dim[i].logarithm && dim[i].normalize)
			{
				for (Dato dato : datos)
				{
					dato.Norm[i] = (dato[i] - dim[i].min) / (dim[i].max - dim[i].min);
				}
			}
				//se logaritma y se normaliza
			else if (dim[i].logarithm && dim[i].normalize)
			{
				for (Dato dato : datos)
				{
					dato.Norm[i] = (Math.log(dato[i]) - Math.log(dim[i].min)) / (Math.log(dim[i].max) - (Math.log(dim[i].min)));
				}
			}
		}
	}

	//endregion

	//region Properties

	/**  Clusters, as a result of the clustering process 
	*/
	public final Cluster[] getClusters()
	{
		return kmc.Clusters;
	}

	//endregion

	//region Clustering

	/**  Main calculation method 
	*/
	public final void Clustering()
	{
		if (catDim != null && catDatos != null && catDim.size() > 0 && catDatos.size() > 0)
		{
			kmc = new KMC(scalDim.toArray(new Dimension[0]), scalDatos.toArray(new Dato[0]), catDim.toArray(new CatDimension[0]), catDatos);
		}
		else
		{
			kmc = new KMC(dim, datos);
		}
		this.dim = kmc.Dim;
		this.datos = kmc.Datos;
		Normalize();
		kmc.Clustering();
		SetResults();
	}

	//endregion

	//region API Obi 

	public final boolean LoadScalDim(Dimension sDim)
	{
		scalDim.add(sDim);
		return true;
	}

	public final boolean LoadCatDim(CatDimension cDim)
	{
		catDim.add(cDim);
		return true;
	}

	public final void LoadData(Map<Integer, ArrayList<Object>> dataPoints)
	{
		HashMap<Integer, HashMap<String, Integer>> catValues = new HashMap<Integer, HashMap<String, Integer>>();
		ArrayList<Object> dataPoint;
		scalDatos = new ArrayList<Dato>();
		catDatos = new ArrayList<ArrayList<boolean[]>>();
		scalDatos = new ArrayList<Dato>();
		Dato sd;
		double[] scalValues;
		HashMap<String, Integer> catValuesForACat;
		for (int i = 0; i < catDim.size(); i++)
		{
			catValues.put(i, new HashMap<String, Integer>());
		}

		for (int index : dataPoints.keySet())
		{
			dataPoint = dataPoints.get(index);

			//variables escalares
			scalValues = new double[scalDim.size()];
			for (int i = 0; i < scalDim.size(); i++)
			{
				scalValues[i] = (Double)dataPoint.get(i);
			}
			sd = new Dato(String.valueOf(index), scalValues);
			scalDatos.add(sd);

			//variables categoricas: primera pasada, categorias
			for (int i = scalDim.size(); i < scalDim.size() + catDim.size(); i++)
			{
				catValuesForACat = catValues.get(i - scalDim.size());
				if (!catValuesForACat.containsKey(dataPoint.get(i).toString()))
				{
					catValuesForACat.put(dataPoint.get(i).toString(), catValuesForACat.size());
				}
			}
		}

		//variables categoricas: segunda pasada, seteo de las categorias
		for (int index : catValues.keySet())
		{
			catDim.get(index).nCats = catValues.get(index).size();
			String[] catNames = new String[catValues.get(index).size()];
			int i = 0;
			for (String cat : catValues.get(index).keySet())
			{
				catNames[i] = cat;
				i++;
			}
			catDim.get(index).catNames = catNames;
			catDim.get(index).catValues = new boolean[catNames.length];
		}

		//variables categoricas: tercera pasada datos de las categorias
		ArrayList<boolean[]> catDato;
		for (int index : dataPoints.keySet())
		{
			dataPoint = dataPoints.get(index);

			catDato = new ArrayList<boolean[]>();
			catDatos.add(catDato);
			for (int i = scalDim.size(); i < scalDim.size() + catDim.size(); i++)
			{
				catDato.add(new boolean[catDim.get(i - scalDim.size()).nCats]);
				for (int j = 0; j < catDato.get(i - scalDim.size()).length; j++)
				{
					catDato.get(i - scalDim.size())[j] = (dataPoint.get(i).toString().equals(catDim.get(i - scalDim.size()).catNames[j]));
				}
			}
		}

		//carga de la api interna
		LoadData(scalDim, scalDatos, catDim, catDatos);
	}

	//catValuesForACat = catValues[i - scalDim.Count];
	//catDim[i - scalDim.Count].SetCategory((string)data[i]);

	//endregion

	//region Log Generation

	private void SetResults()
	{
		//sort
		for (Cluster clust : getClusters())
		{
			clust.Indice = dim[0].weight * clust.Centro.Valores[0] * (dim[0].sort? 1 : 0);
			if (dim.length > 1)
			{
				clust.Indice += dim[1].weight * clust.Centro.Valores[1] * (dim[1].sort? 1 : 0);
			}
			else if (dim.length > 2)
			{
				clust.Indice += dim[2].weight * clust.Centro.Valores[2] * (dim[2].sort? 1 : 0);
			}
		}
		resultClusters = new ArrayList(kmc.Clusters);
		Collections.sort(resultClusters);

		int indice = 1;
		for (Cluster clust : getClusters())
		{
			clust.Indice = indice++;
		}
	}

	/**  Generate a speciphic log for ECC 
	*/
	public final void GenerateECCLog()
	{
		SetResults();
		try
		{
//C# TO JAVA CONVERTER WARNING: The java.io.OutputStreamWriter constructor does not accept all the arguments passed to the System.IO.StreamWriter constructor:
//ORIGINAL LINE: StreamWriter ecc = new StreamWriter("..\\ECC_KMC.txt",false);
			java.io.OutputStreamWriter ecc = new java.io.OutputStreamWriter("..\\ECC_KMC.txt");
			double dist;
			for (Dato dato : this.datos)
			{
				ecc.write(String.valueOf(dato.Codigo.PadRight(15,' ')));
				ecc.write(String.valueOf(tangible.DotNetToJavaStringHelper.padLeft(dato.Cluster.Indice.toString(), 4,' ')));
				dist = dato[0] * dim[0].weight + dato[1] * dim[1].weight;
				ecc.write(String.valueOf(tangible.DotNetToJavaStringHelper.padLeft(String.valueOf(dist), 18,' ')));
				ecc.write(";" + System.lineSeparator());
			}
			ecc.close();
			System.Diagnostics.Debug.WriteLine("Completado archivo KMC.txt");
		}
		catch (java.lang.Exception e)
		{
			System.Diagnostics.Debug.WriteLine("Error");
		}
	}

	/**  Generate log 
	*/
	public final void GenerateLog()
	{
		SetResults();
		try
		{
//C# TO JAVA CONVERTER WARNING: The java.io.OutputStreamWriter constructor does not accept all the arguments passed to the System.IO.StreamWriter constructor:
//ORIGINAL LINE: StreamWriter so = new StreamWriter("..\\KMC.txt",false);
			java.io.OutputStreamWriter so = new java.io.OutputStreamWriter("..\\KMC.txt");
			so.write("CLUSTERS:" + System.lineSeparator());
			so.write("---------" + System.lineSeparator());
			so.write("" + System.lineSeparator());
			for (Cluster clust : resultClusters)
			{
				so.write(clust.Indice + "\t" + clust.Centro.Valores[0].toString("//0.00") + "\t" + clust.Centro.Valores[1].toString("//0.00") + "\t" + clust.Centro.Valores[2].toString("//0.00") + System.lineSeparator());
			}
			so.write("\nDATOS:" + System.lineSeparator());
			so.write("------" + System.lineSeparator());
			so.write("" + System.lineSeparator());
			for (Cluster clust : resultClusters)
			{
				 for (Dato dato : clust.Datos)
				 {
					 so.write(String.valueOf(dato.Codigo + "\t" + dato[0].toString("//0.00") + "\t" + dato[1].toString("//0.00") + "\t" + dato[2].toString("//0.00") + "\t" + clust.Indice) + System.lineSeparator());
				 }
			}
			so.close();
			System.Diagnostics.Debug.WriteLine("Completado archivo KMC.txt");
		}
		catch (java.lang.Exception e)
		{
			System.Diagnostics.Debug.WriteLine("Error");
		}
	}

	//endregion

}