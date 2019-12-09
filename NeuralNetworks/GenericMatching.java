package NeuralNetworks;

import java.util.*;

//region Imports


//endregion


/**  Class for matching Generics to components/ 
*/
public class GenericMatching
{

	//region Fields

	private String generic;
	private Red matchingNN;
	private HashMap<String, ArrayList<String>> rules;
	private ArrayList<String> components;
	private HashMap<String, Integer> inputs;
	private HashMap<String, Output> outputs;

	private double m;
	private double n;
	private double interProp;

	//endregion

	//region Constructor

	/**  Constructor 
	 @param generic the generic code 
	 @param m learning coefficient 
	 @param n momentum factor 
	 @param interProp proportion of neurons in intermediate layer with respect of input layer 
	*/
	public GenericMatching(String generic, double m, double n, double interProp)
	{
		this.generic = generic;
		this.m = m;
		this.n = n;
		this.interProp = interProp;
	}

	//endregion

	//region Public Methods

	/**  Load data for calculation 
	 @param rules features and their options 
	 @param components list of components for this generic 
	*/
	public final void LoadData(HashMap<String, ArrayList<String>> rules, ArrayList<String> components)
	{

		this.inputs = new HashMap<String, Integer>();
		int i = 0;
		for (String feature : rules.keySet())
		{
			for (int j = 0; j < rules.get(feature).size(); j++)
			{
				this.inputs.put(feature + "|" + rules.get(feature).get(j), i);
				i++;
			}
		}

		outputs = new HashMap<String,Output>();
		Output output;
		String code;
		for (int j = 0; j < components.size(); j++)
		{
			code = components.get(j);
			output = new Output(code);
			output.index = j;
			outputs.put(code, output);
		}
		int interCount = (int)(inputs.size() * interProp);
		matchingNN = new Red(1, inputs.size(), interCount, outputs.size(), 0.1, 0.2, Neurona.Funcion.sigmoidal);
	}

	/**  Calculate selected components for selected options 
	 @param selOptions selected options for each feature 
	 @return  selected components and their quantities 
	*/
	public final HashMap<String, Double> Calculate(HashMap<String, String> selOptions)
	{
		HashMap<String, Double> selComponents = new HashMap<String, Double>();
		double[] inpt = InputToArray(selOptions);
		matchingNN.procesar(inpt);
		double[] outt = matchingNN.getSalida();
		return ArrayToOutput(outt);
	}

	/**  Train the system 
	 @param selOptions selected options for each feacutre 
	 @param selComponents selected components and their quantities 
	*/
	public final void Train(HashMap<String, String> selOptions, HashMap<String, Double> selComponents)
	{
		for (String compCode : selComponents.keySet())
		{
			if (outputs.containsKey(compCode) && selComponents.get(compCode).compareTo(outputs.get(compCode).max) > 0)
			{
				outputs.get(compCode).max = selComponents.get(compCode);
			}
		}
		double[] inpt = InputToArray(selOptions);
		double[] outt = OutputToArray(selComponents);
		matchingNN.entrenar(inpt, outt);
	}

	/**  Save neural network weights 
	*/
	public final void Save()
	{
//C# TO JAVA CONVERTER TODO TASK: C# to Java Converter cannot determine whether this System.IO.FileStream is input or output:
		FileStream weights = new FileStream("..\\..\\..\\weights_" + generic + ".txt", FileMode.Create);
		SaveTo(weights);
	}

	/**  Read neural network weights 
	*/
	public final void Read()
	{
//C# TO JAVA CONVERTER TODO TASK: C# to Java Converter cannot determine whether this System.IO.FileStream is input or output:
		FileStream weights = new FileStream("..\\..\\..\\weights_" + generic + ".txt", FileMode.Create);
		ReadFrom(weights);
	}

	/**  Save to a file stream 
	 @param weights stream of weights 
	*/
	public final void SaveTo(java.io.OutputStream weights)
	{
		String linea;
		Neurona ne;
		int indice;
		for (int i = 0; i < matchingNN.Oculta; i++)
		{
			ne = matchingNN.getNeurona(1, i);
			linea = 1 + " " + i + " " + matchingNN.Entrada + "";
			for (indice = 0; indice < linea.length(); indice++)
			{
//C# TO JAVA CONVERTER WARNING: Unsigned integer types have no direct equivalent in Java:
//ORIGINAL LINE: weights.WriteByte((byte)linea[indice]);
				weights.write((byte)linea.charAt(indice));
			}
//C# TO JAVA CONVERTER WARNING: Unsigned integer types have no direct equivalent in Java:
//ORIGINAL LINE: weights.WriteByte((byte)' ');
			weights.write((byte)' ');
			for (int j = 0; j < matchingNN.Entrada; j++)
			{
				linea = ne.getW(j) + "";
				for (indice = 0; indice < linea.length(); indice++)
				{
//C# TO JAVA CONVERTER WARNING: Unsigned integer types have no direct equivalent in Java:
//ORIGINAL LINE: weights.WriteByte((byte)linea[indice]);
					weights.write((byte)linea.charAt(indice));
				}
//C# TO JAVA CONVERTER WARNING: Unsigned integer types have no direct equivalent in Java:
//ORIGINAL LINE: weights.WriteByte((byte)' ');
				weights.write((byte)' ');
			}
//C# TO JAVA CONVERTER WARNING: Unsigned integer types have no direct equivalent in Java:
//ORIGINAL LINE: weights.WriteByte((byte)'\r');
			weights.write((byte)'\r');
//C# TO JAVA CONVERTER WARNING: Unsigned integer types have no direct equivalent in Java:
//ORIGINAL LINE: weights.WriteByte((byte)'\n');
			weights.write((byte)'\n');
		}

		for (int i = 0; i < matchingNN.Salida; i++)
		{
			ne = matchingNN.getNeurona(2, i);
			linea = "2" + " " + i + " " + matchingNN.Oculta + " ";
			for (indice = 0; indice < linea.length(); indice++)
			{
//C# TO JAVA CONVERTER WARNING: Unsigned integer types have no direct equivalent in Java:
//ORIGINAL LINE: weights.WriteByte((byte)linea[indice]);
				weights.write((byte)linea.charAt(indice));
			}
			for (int j = 0; j < matchingNN.Oculta; j++)
			{
				linea = ne.getW(j) + " ";
				for (indice = 0; indice < linea.length(); indice++)
				{
//C# TO JAVA CONVERTER WARNING: Unsigned integer types have no direct equivalent in Java:
//ORIGINAL LINE: weights.WriteByte((byte)linea[indice]);
					weights.write((byte)linea.charAt(indice));
				}
			}
		}
//C# TO JAVA CONVERTER WARNING: Unsigned integer types have no direct equivalent in Java:
//ORIGINAL LINE: weights.WriteByte((byte)'f');
		weights.write((byte)'f');
//C# TO JAVA CONVERTER WARNING: Unsigned integer types have no direct equivalent in Java:
//ORIGINAL LINE: weights.WriteByte((byte)'\r');
		weights.write((byte)'\r');
//C# TO JAVA CONVERTER WARNING: Unsigned integer types have no direct equivalent in Java:
//ORIGINAL LINE: weights.WriteByte((byte)'\n');
		weights.write((byte)'\n');

	}

	/**  Read from a file stream 
	 @param weights stream of weights 
	*/
	public final void ReadFrom(java.io.InputStream weights)
	{
		char caracter;
		ArrayList<String> lineas = new ArrayList<String>();
		String linea = "";
		int byteInt;
//C# TO JAVA CONVERTER WARNING: Unsigned integer types have no direct equivalent in Java:
//ORIGINAL LINE: int fin = Convert.ToByte('f');
		int fin = (byte)'f';
		while (true)
		{
			byteInt = weights.read();
			if (byteInt == fin)
			{
				//Debug.WriteLine("[" + linea + "]");
				lineas.add(linea);
				break;
			}
			caracter = (char)byteInt;
			if (caracter == '\r')
			{
				weights.read();
				//Debug.WriteLine("[" + linea + "]");
				lineas.add(linea);
				linea = "";
			}
			else
			{
				linea += caracter;
			}
		}
		weights.close();

		char[] separador = {' '};
		String[] valores;
		for (String lin : lineas)
		{
			valores = lin.split(java.util.regex.Pattern.quote(separador.toString()), -1);
			int capa = Integer.parseInt(valores[0]);
			int indice = Integer.parseInt(valores[1]);
			int cantConexiones = Integer.parseInt(valores[2]);
			for (int i = 0; i < cantConexiones; i++)
			{
				matchingNN.getNeurona(capa, indice).setW(Double.parseDouble(valores[i + 3]), i);
			}
		}
	}

	//endregion

	//region Private Methods

	//region Auxiliar Methods

	private double[] InputToArray(HashMap<String, String> selOptions)
	{
		double[] inpt = new double[inputs.size()];
		for (String feature : selOptions.keySet())
		{
			String option = selOptions.get(feature);
			int index = inputs.get(feature + "|" + option);
			inpt[index] = 1.0;
		}
		return inpt;
	}

	private double[] OutputToArray(HashMap<String, Double> selComponents)
	{
		double[] outt = new double[outputs.size()];
		for (String selComp : selComponents.keySet())
		{
			if (!outputs.containsKey(selComp))
			{
				throw new RuntimeException("Error. Component " + selComp + " does not exist.");
			}
			int index = outputs.get(selComp).index;
			outt[index] = 0;
			if (outputs.get(selComp).max > 0)
			{
				outt[index] = selComponents.get(selComp) / outputs.get(selComp).max;
			}

		}
		return outt;
	}

	private HashMap<String, Double> ArrayToOutput(double[] outt)
	{
		HashMap<String, Double> output = new HashMap<String, Double>();

		double val;
		for (Output outp : outputs.values())
		{
			val = Math.round(outt[outp.index] * outp.max);
			if (val > 0)
			{
				output.put(outp.code, val);
			}
		}
		return output;
	}

	//endregion

	//endregion

	//region Class output

	public static class Output
	{
		public String code;
		public int index;
		public double max;

		public Output(String code)
		{
			this.code = code;
			this.index = 0;
			this.max = 0;
		}
	}

	//endregion
}