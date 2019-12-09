package Clustering;

import ClassicalStat.*;
import NumCalc.*;
import java.util.*;

//region Imports


//endregion


public class SkuClustering
{

	//region Fields

	private Statistics stat;
	private Functions func;
	private PCA pca;
	private ECCClustering clust;

	private HashMap<Integer, Color> colorsByIndex;
	private HashMap<String, Color> colorsByName;
	private ArrayList<Sku> skus;

	private int nDim;
	private int nData;
	private double[][] data;
	private double[][] s;
	private double[] v;
	private double[][] lts;
	private HashMap<String, PropertyInfo> propsDict;
	private ArrayList<PropertyInfo> props;
	private double percInfo;
	private int nSelDims;

	private Dimension[] scalDimensions;
	private Dato[] scalData;
	private ArrayList<Cluster> clusters;
	private HashMap<String, Dato> scalDataDict;
	private HashMap<String, Sku> skusDict;


	//endregion

	//region Constructor

	public SkuClustering()
	{
		stat = new Statistics();
		func = new Functions();
		pca = new PCA();
		clust = new ECCClustering();

		colorsByIndex = new HashMap<Integer, Color>();
		colorsByName = new HashMap<String, Color>();
		LoadColors();

		propsDict = new HashMap<String, PropertyInfo>();
		props = new ArrayList<PropertyInfo>();
		scalDataDict = new HashMap<String, Dato>();
		skusDict = new HashMap<String, Sku>();

		nDim = 10;
		nData = 100;
		data = new double[nData][nDim];
		percInfo = 90;
		nSelDims = 1;
	}

	//endregion

	//region Properties

	public final double[][] getData()
	{
		return data;
	}

	public final double[][] getScores()
	{
		return pca.Scores;
	}

	public final double[] getVars()
	{
		return pca.Vars;
	}

	public final double getPercInfo()
	{
		return percInfo;
	}
	public final void setPercInfo(double value)
	{
		percInfo = value;
	}

	public final int getNSelDims()
	{
		return nSelDims;
	}
	public final void setNSelDims(int value)
	{
		nSelDims = value;
	}

	public final ArrayList<Cluster> getClusters()
	{
		return clusters;
	}

	//endregion

	//region Public Methods

	//region Main Methods

	public final void LoadData(ArrayList<Sku> skus)
	{
		this.skus = skus;
		for (Sku sku : skus)
		{
			skusDict.put(sku.getcode(), sku);
		}
		this.nData = skus.size();

		java.lang.Class type = Sku.class;
		for (PropertyInfo pi : type.GetProperties())
		{
			if (pi.PropertyType != Double.class)
			{
				continue;
			}
			propsDict.put(pi.Name, pi);
			props.add(pi);
		}
		this.nDim = propsDict.size();

		data = new double[nData][nDim];
		int i;
		for (int j = 0; j < skus.size(); j++)
		{
			i = 0;
			for (PropertyInfo pi : propsDict.values())
			{
				data[j][i++] = (Double)pi.GetValue(skus.get(j), null);
			}
		}
	}

	public final HashMap<String, Double> CalcCorrelations(String demand)
	{
		HashMap<String, Double> propCorrelations = new HashMap<String,Double>();
		java.lang.Class type = Sku.class;
		ArrayList<Double> demVals = new ArrayList<Double>();
		PropertyInfo dem = type.GetProperty(demand);
		for (Sku sku : skus)
		{
			demVals.add((Double)dem.GetValue(sku, null));
		}

		ArrayList<ArrayList<Double>> propVals = new ArrayList<ArrayList<Double>>();
		int n = 0;
		double r2;
		for (PropertyInfo pi : propsDict.values())
		{
			propVals.add(new ArrayList<Double>());
			for (int j = 0;j < skus.size();j++)
			{
				propVals.get(n).add((Double)pi.GetValue(skus.get(j), null));
			}
			r2 = stat.R2(propVals.get(n), demVals);
			propCorrelations.put(pi.Name, r2);
			n++;
		}

		return propCorrelations;
	}

	public final void CalcPCA()
	{
		pca.Calculate(data, nData, nDim);
		double totVars = 0;
		double acumVars = 0;
		for (int i = 0; i < pca.Vars.getLength();i++)
		{
			totVars += pca.Vars[i];
		}
		for (int i = 0; i < pca.Vars.getLength(); i++)
		{
			acumVars += pca.Vars[i] / totVars;
			if (acumVars > percInfo)
			{
				nSelDims = i - 1;
				return;
			}
		}
	}

	public final double[][] GetLTValues(int nSelDims)
	{
		lts = pca.GetLinearTransfs(nSelDims);
		return lts;
	}

	public final void Clustering()
	{
		lts = pca.GetLinearTransfs(nSelDims);
		clust = new ECCClustering();
		LoadScalDim();
		LoadScalData();
		clust.LoadData(scalDimensions, scalData);
		clust.Clustering();
		clusters = new ArrayList<Cluster>(clust.Clusters);
	}

	public final HashMap<String, Sku> SelectSimilSkus(Sku sku, int nSimils)
	{
		HashMap<String, Sku> simils = new HashMap<String, Sku>();
		Dato dato = GetDato(sku, nSelDims);
		double dist;
		double minDist = Double.MAX_VALUE;
		Cluster skuCluster = null;
		for (Cluster c : clusters)
		{
			dist = func.DistEuclid(c.Centro.Valores, dato.Valores);
			if (dist < minDist)
			{
				minDist = dist;
				skuCluster = c;
			}
		}
		ArrayList<Dato> skuData = new ArrayList<Dato>();
		for (Dato d : skuCluster.Datos)
		{
			d.Dist = func.DistEuclid(d.Valores, dato.Valores);
			skuData.add(d);
		}
		Collections.sort(skuData, Dato.sortDistAsc());
		double totWeight = 0;
		for (int i = 0; i < nSimils; i++)
		{
			Sku skU = skusDict.get(skuData.get(i).Codigo);
			if (skuData.get(i).Dist == 0)
			{
				skU.setsimilWeight(1000000);
			}
			else
			{
				skU.setsimilWeight(1.0 / skuData.get(i).Dist);
			}
			totWeight += skU.getsimilWeight();
			simils.put(skU.getcode(), skU);
		}
		for (Sku skU : simils.values())
		{
			skU.setsimilWeight(skU.getsimilWeight() / totWeight);
		}

		return simils;
	}

	//endregion

	//region Auxiliar Methods

	public final Color GetColor(int index)
	{
		if (!colorsByIndex.containsKey(index))
		{
			throw new RuntimeException("Error. Color not defined");
		}
		return colorsByIndex.get(index);
	}

	public final Color GetColor(String name)
	{
		if (!colorsByName.containsKey(name))
		{
			throw new RuntimeException("Error. Color not defined");
		}
		return colorsByName.get(name);
	}

	//endregion

	//endregion

	//region Private Methods

	private double[] GetPropArray(Sku sku)
	{
		double[] propArray = new double[propsDict.size()];
		for (int i = 0; i < propsDict.size(); i++)
		{
			propArray[i] = (Double)props.get(i).GetValue(sku, null);
		}
		return propArray;
	}

	private Dato GetDato(Sku sku, int nSelDims)
	{
		double[] propArray = GetPropArray(sku);
		double[] lt = pca.GetLinearTransf(propArray, nSelDims);
		Dato dato = new Dato(sku.getcode(), lt);
		return dato;
	}

	private void LoadScalDim()
	{
		scalDimensions = new Dimension[nSelDims];
		Dimension dim;
		double[] min = new double[nSelDims];
		for (int i = 0; i < min.length; i++)
		{
			min[i] = Double.MAX_VALUE;
		}
		double[] max = new double[nSelDims];
		for (int i = 0; i < max.length; i++)
		{
			max[i] = -Double.MAX_VALUE;
		}
		for (int i = 0; i < skus.size(); i++)
		{
			for (int j = 0; j < nSelDims; j++)
			{
				if (lts[i][j] < min[j])
				{
					min[j] = lts[i][j];
				}
				if (lts[i][j] > max[j])
				{
					max[j] = lts[i][j];
				}
			}
		}
		for (int i = 0; i < nSelDims; i++)
		{
			dim = new Dimension();
			dim.index = i;
			dim.min = min[i];
			dim.max = max[i];
			dim.logarithm = false;
			dim.normalize = true;
			dim.nDivisions = 5;
			dim.sort = false;
			dim.weight = 1;
			scalDimensions[i] = dim;
		}
	}

	private void LoadScalData()
	{
		this.scalDataDict = new HashMap<String, Dato>();
		Dato dato;
		double[] values;
		scalData = new Dato[skus.size()];
		for (int i = 0; i < skus.size(); i++)
		{
			values = new double[nSelDims];
			for (int j = 0; j < nSelDims; j++)
			{
				values[j] = lts[i][j];
			}
			dato = new Dato(skus.get(i).getcode(), values);
			scalData[i] = dato;
			scalDataDict.put(dato.Codigo, dato);
		}
	}

	private Dato GetScalDatum(Sku sku)
	{
		Dato dato;
		double[] values;
		scalData = new Dato[skus.size()];
		values = new double[nSelDims + 1];
		for (int i = 0; i <= nSelDims; i++)
		{
		   //values[i] = sku.pca.Scores[i, j];
		}
		dato = new Dato(sku.getcode(), values);
		return dato;
	}

	private void LoadColors()
	{
		Color c;
		c = new Color(1, "IndianRed", ColorGroup.Red, 205, 92, 92);
		AddColor(c);
		c = new Color(2, "LightCoral", ColorGroup.Red, 240, 128, 128);
		AddColor(c);
		c = new Color(3, "Salmon", ColorGroup.Red, 250, 128, 114);
		AddColor(c);
		c = new Color(4, "DarkSalmon", ColorGroup.Red, 233, 150, 122);
		AddColor(c);
		c = new Color(5, "Crimson", ColorGroup.Red, 220, 20, 60);
		AddColor(c);
		c = new Color(6, "Red", ColorGroup.Red, 265, 0, 0);
		AddColor(c);
		c = new Color(7, "FireBrick", ColorGroup.Red, 178, 34, 34);
		AddColor(c);
		c = new Color(8, "DarkRed", ColorGroup.Red, 139, 0, 0);
		AddColor(c);
		c = new Color(9, "Pink", ColorGroup.Rose, 255, 192, 203);
		AddColor(c);
		c = new Color(10, "LightPink", ColorGroup.Rose, 255, 182, 193);
		AddColor(c);
		c = new Color(11, "HotPink", ColorGroup.Rose, 255, 105, 180);
		AddColor(c);
		c = new Color(12, "DeepPink", ColorGroup.Rose, 255, 20, 147);
		AddColor(c);
		c = new Color(13, "MediumVioletRed", ColorGroup.Rose, 199, 21, 133);
		AddColor(c);
		c = new Color(14, "PaleVioletRed", ColorGroup.Rose, 219, 112, 147);
		AddColor(c);
		c = new Color(15, "LightSalmon", ColorGroup.Orange, 255, 160, 122);
		AddColor(c);
		c = new Color(16, "Coral", ColorGroup.Orange, 255, 127, 80);
		AddColor(c);
		c = new Color(17, "Tomato", ColorGroup.Orange, 255, 99, 71);
		AddColor(c);
		c = new Color(18, "OrangeRed", ColorGroup.Orange, 255, 69, 0);
		AddColor(c);
		c = new Color(19, "DarkOrange", ColorGroup.Orange, 255, 140, 0);
		AddColor(c);
		c = new Color(20, "Orange", ColorGroup.Orange, 255, 165, 0);
		AddColor(c);
		c = new Color(21, "Gold", ColorGroup.Yellow, 255, 215, 0);
		AddColor(c);
		c = new Color(22, "Yellow", ColorGroup.Yellow, 255, 255, 0);
		AddColor(c);
		c = new Color(23, "LightYellow", ColorGroup.Yellow, 255, 255, 224);
		AddColor(c);
		c = new Color(24, "LemmonChiffon", ColorGroup.Yellow, 255, 250, 205);
		AddColor(c);
		c = new Color(25, "LightGoldenrodYellow", ColorGroup.Yellow, 250, 250, 210);
		AddColor(c);
		c = new Color(26, "PapayaWhip", ColorGroup.Yellow, 255, 239, 213);
		AddColor(c);
		c = new Color(27, "Moccasin", ColorGroup.Yellow, 255, 228, 181);
		AddColor(c);
		c = new Color(28, "PeachPuff", ColorGroup.Yellow, 2155, 218, 185);
		AddColor(c);
		c = new Color(29, "PaleGoldenrod", ColorGroup.Yellow, 238, 232, 170);
		AddColor(c);
		c = new Color(30, "Khaki", ColorGroup.Yellow, 240, 230, 140);
		AddColor(c);
		c = new Color(31, "DarkKhaki", ColorGroup.Yellow, 189, 183, 107);
		AddColor(c);
		c = new Color(32, "Lavender", ColorGroup.Purple, 230, 230, 250);
		AddColor(c);
		c = new Color(33, "Thistle", ColorGroup.Purple, 216, 191, 216);
		AddColor(c);
		c = new Color(34, "Plum", ColorGroup.Purple, 221, 160, 221);
		AddColor(c);
		c = new Color(35, "Violet", ColorGroup.Purple, 238, 130, 238);
		AddColor(c);
		c = new Color(36, "Orchid", ColorGroup.Purple, 218, 112, 214);
		AddColor(c);
		c = new Color(37, "Fucsia", ColorGroup.Purple, 255, 0, 255);
		AddColor(c);
		c = new Color(38, "MediumOrchid", ColorGroup.Purple, 186, 85, 211);
		AddColor(c);
		c = new Color(39, "MediumPurple", ColorGroup.Purple, 147, 112, 219);
		AddColor(c);
		c = new Color(40, "BlueViolet", ColorGroup.Purple, 138, 43, 226);
		AddColor(c);
		c = new Color(41, "DarkViolet", ColorGroup.Purple, 148, 0, 211);
		AddColor(c);
		c = new Color(42, "DarkOrchid", ColorGroup.Purple, 153, 50, 204);
		AddColor(c);
		c = new Color(43, "DarkMagenta", ColorGroup.Purple, 139, 0, 139);
		AddColor(c);
		c = new Color(44, "Purple", ColorGroup.Purple, 128, 0, 128);
		AddColor(c);
		c = new Color(45, "Indigo", ColorGroup.Purple, 75, 0, 130);
		AddColor(c);
		c = new Color(46, "StateBlue", ColorGroup.Purple, 106, 90, 205);
		AddColor(c);
		c = new Color(47, "DarkStateBlue", ColorGroup.Purple, 72, 61, 139);
		AddColor(c);
		c = new Color(48, "GreenYellow", ColorGroup.Green, 173, 255, 47);
		AddColor(c);
		c = new Color(49, "Chartreuse", ColorGroup.Green, 127, 255, 0);
		AddColor(c);
		c = new Color(50, "LawnGreen", ColorGroup.Green, 124, 252, 0);
		AddColor(c);
		c = new Color(51, "Lime", ColorGroup.Green, 0, 255, 0);
		AddColor(c);
		c = new Color(52, "LimeGreen", ColorGroup.Green, 50, 205, 50);
		AddColor(c);
		c = new Color(53, "PaleGreen", ColorGroup.Green, 152, 251, 152);
		AddColor(c);
		c = new Color(54, "LightGreen", ColorGroup.Green, 144, 238, 144);
		AddColor(c);
		c = new Color(55, "MediumSpringGreen", ColorGroup.Green, 0, 250, 154);
		AddColor(c);
		c = new Color(56, "SpringGreen", ColorGroup.Green, 0, 255, 127);
		AddColor(c);
		c = new Color(57, "MediumSeaGreen", ColorGroup.Green, 60, 179, 113);
		AddColor(c);
		c = new Color(58, "SeaGreen", ColorGroup.Green, 46, 139, 87);
		AddColor(c);
		c = new Color(59, "ForestGreen", ColorGroup.Green, 34, 139, 34);
		AddColor(c);
		c = new Color(60, "Green", ColorGroup.Green, 0, 128, 0);
		AddColor(c);
		c = new Color(61, "DarkGreen", ColorGroup.Green, 0, 100, 0);
		AddColor(c);
		c = new Color(62, "YellowGreen", ColorGroup.Green, 154, 205, 50);
		AddColor(c);
		c = new Color(63, "OliveDrav", ColorGroup.Green, 107, 142, 35);
		AddColor(c);
		c = new Color(64, "Olive", ColorGroup.Green, 128, 128, 0);
		AddColor(c);
		c = new Color(65, "DarkOliveGreen", ColorGroup.Green, 85, 107, 47);
		AddColor(c);
		c = new Color(66, "MediumAquamarine", ColorGroup.Green, 102, 205, 170);
		AddColor(c);
		c = new Color(67, "DarkSeaGreen", ColorGroup.Green, 143, 188, 143);
		AddColor(c);
		c = new Color(68, "LightSeaGreen", ColorGroup.Green, 32, 178, 170);
		AddColor(c);
		c = new Color(69, "DarkCyan", ColorGroup.Green, 0, 139, 139);
		AddColor(c);
		c = new Color(70, "Teal", ColorGroup.Green, 0, 128, 128);
		AddColor(c);
		c = new Color(71, "AquaCyan", ColorGroup.Blue, 0, 255, 255);
		AddColor(c);
		c = new Color(72, "LightCyan", ColorGroup.Blue, 224, 255, 255);
		AddColor(c);
		c = new Color(73, "PaleTurquoise", ColorGroup.Blue, 175, 238, 238);
		AddColor(c);
		c = new Color(74, "Aquamarine", ColorGroup.Blue, 127, 255, 212);
		AddColor(c);
		c = new Color(75, "Turquoise", ColorGroup.Blue, 64, 224, 208);
		AddColor(c);
		c = new Color(76, "MediumTurquoise", ColorGroup.Blue, 72, 209, 204);
		AddColor(c);
		c = new Color(77, "DarkTurquoise", ColorGroup.Blue, 0, 206, 209);
		AddColor(c);
		c = new Color(78, "CadetBlue", ColorGroup.Blue, 95, 158, 160);
		AddColor(c);
		c = new Color(79, "SteelBlue", ColorGroup.Blue, 70, 130, 180);
		AddColor(c);
		c = new Color(80, "LightSteelBlue", ColorGroup.Blue, 176, 196, 222);
		AddColor(c);
		c = new Color(81, "PowderBlue", ColorGroup.Blue, 176, 224, 230);
		AddColor(c);
		c = new Color(82, "LightBlue", ColorGroup.Blue, 173, 216, 230);
		AddColor(c);
		c = new Color(83, "SkyBlue", ColorGroup.Blue, 135, 206, 235);
		AddColor(c);
		c = new Color(84, "LightSkyBlue", ColorGroup.Blue, 135, 206, 250);
		AddColor(c);
		c = new Color(85, "DeepSkyBlue", ColorGroup.Blue, 0, 191, 255);
		AddColor(c);
		c = new Color(86, "DodgerBlue", ColorGroup.Blue, 30, 144, 255);
		AddColor(c);
		c = new Color(87, "CornFlowerBlue", ColorGroup.Blue, 100, 149, 237);
		AddColor(c);
		c = new Color(88, "MediumSlateBlue", ColorGroup.Blue, 123, 104, 238);
		AddColor(c);
		c = new Color(89, "RoyalBlue", ColorGroup.Blue, 65, 105, 225);
		AddColor(c);
		c = new Color(90, "Blue", ColorGroup.Blue, 0, 0, 255);
		AddColor(c);
		c = new Color(91, "MediumBlue", ColorGroup.Blue, 0, 0, 205);
		AddColor(c);
		c = new Color(92, "DarkBlue", ColorGroup.Blue, 0, 0, 139);
		AddColor(c);
		c = new Color(93, "Navy", ColorGroup.Blue, 0, 0, 128);
		AddColor(c);
		c = new Color(94, "MidnightBlue", ColorGroup.Blue, 25, 25, 112);
		AddColor(c);
		c = new Color(95, "ComSick", ColorGroup.Brown, 255, 248, 220);
		AddColor(c);
		c = new Color(96, "BlanchedAlmond", ColorGroup.Brown, 255, 235, 205);
		AddColor(c);
		c = new Color(97, "Bisque", ColorGroup.Brown, 255, 228, 196);
		AddColor(c);
		c = new Color(98, "NavajoWhite", ColorGroup.Brown, 255, 222, 173);
		AddColor(c);
		c = new Color(99, "Wheat", ColorGroup.Brown, 245, 222, 179);
		AddColor(c);
		c = new Color(100, "BurlyWood", ColorGroup.Brown, 222, 184, 135);
		AddColor(c);
		c = new Color(101, "Tan", ColorGroup.Brown, 210, 180, 140);
		AddColor(c);
		c = new Color(102, "RosyBrown", ColorGroup.Brown, 188, 143, 143);
		AddColor(c);
		c = new Color(103, "SandyBrown", ColorGroup.Brown, 244, 164, 96);
		AddColor(c);
		c = new Color(104, "Goldenrod", ColorGroup.Brown, 218, 165, 32);
		AddColor(c);
		c = new Color(105, "DarkGoldenrod", ColorGroup.Brown, 184, 134, 11);
		AddColor(c);
		c = new Color(106, "Peru", ColorGroup.Brown, 205, 133, 63);
		AddColor(c);
		c = new Color(107, "Chocolate", ColorGroup.Brown, 210, 105, 30);
		AddColor(c);
		c = new Color(108, "SaddleBrown", ColorGroup.Brown, 139, 69, 19);
		AddColor(c);
		c = new Color(109, "Siena", ColorGroup.Brown, 160, 82, 45);
		AddColor(c);
		c = new Color(110, "Brown", ColorGroup.Brown, 165, 42, 42);
		AddColor(c);
		c = new Color(111, "Maroon", ColorGroup.Brown, 128, 0, 0);
		AddColor(c);
		c = new Color(112, "White", ColorGroup.White, 255, 255, 255);
		AddColor(c);
		c = new Color(113, "Snow", ColorGroup.White, 255, 250, 250);
		AddColor(c);
		c = new Color(114, "HoneyDew", ColorGroup.White, 240, 255, 240);
		AddColor(c);
		c = new Color(115, "MintCream", ColorGroup.White, 245, 255, 250);
		AddColor(c);
		c = new Color(116, "Azure", ColorGroup.White, 240, 255, 255);
		AddColor(c);
		c = new Color(117, "AliceBlue", ColorGroup.White, 240, 248, 255);
		AddColor(c);
		c = new Color(118, "GhostWhite", ColorGroup.White, 248, 248, 255);
		AddColor(c);
		c = new Color(119, "WhiteSmoke", ColorGroup.White, 245, 245, 245);
		AddColor(c);
		c = new Color(120, "SeaShell", ColorGroup.White, 255, 245, 238);
		AddColor(c);
		c = new Color(121, "Beige", ColorGroup.White, 245, 245, 220);
		AddColor(c);
		c = new Color(122, "OldLace", ColorGroup.White, 253, 245, 230);
		AddColor(c);
		c = new Color(123, "FloralWhite", ColorGroup.White, 255, 250, 240);
		AddColor(c);
		c = new Color(124, "Ivory", ColorGroup.White, 255, 255, 240);
		AddColor(c);
		c = new Color(125, "AntiqueWhite", ColorGroup.White, 250, 235, 215);
		AddColor(c);
		c = new Color(126, "Linen", ColorGroup.White, 250, 240, 230);
		AddColor(c);
		c = new Color(127, "LavenderBlush", ColorGroup.White, 255, 240, 245);
		AddColor(c);
		c = new Color(128, "MistyRose", ColorGroup.White, 255, 228, 225);
		AddColor(c);
		c = new Color(129, "Gainsboro", ColorGroup.Gray, 220, 220, 220);
		AddColor(c);
		c = new Color(130, "LightGray", ColorGroup.Gray, 211, 211, 211);
		AddColor(c);
		c = new Color(131, "Silver", ColorGroup.Gray, 192, 192, 192);
		AddColor(c);
		c = new Color(132, "DarkGray", ColorGroup.Gray, 169, 169, 169);
		AddColor(c);
		c = new Color(133, "Gray", ColorGroup.Gray, 128, 128, 128);
		AddColor(c);
		c = new Color(134, "DimGray", ColorGroup.Gray, 105, 105, 105);
		AddColor(c);
		c = new Color(135, "LightSlateGray", ColorGroup.Gray, 119, 136, 153);
		AddColor(c);
		c = new Color(136, "SlateGray", ColorGroup.Gray, 112, 128, 144);
		AddColor(c);
		c = new Color(137, "DarkSlateGray", ColorGroup.Gray, 49, 79, 79);
		AddColor(c);
		c = new Color(138, "Black", ColorGroup.Gray, 0, 0, 0);
		AddColor(c);
	}

	private void AddColor(Color c)
	{
		colorsByIndex.put(c.index, c);
		colorsByName.put(c.name, c);
	}

	//endregion

	//region Inner Classes

	public static class Color
	{

		public int index;
		public String name;
		public ColorGroup group = ColorGroup.values()[0];
		public int R;
		public int G;
		public int B;

		public Color(int index, String name, ColorGroup group, int R, int G, int B)
		{
			this.index = index;
			this.name = name;
			this.group = group;
			this.R = R;
			this.G = G;
			this.B = B;
		}

	}

	public enum ColorGroup
	{
		Red,
		Rose,
		Orange,
		Yellow,
		Purple,
		Green,
		Blue,
		Brown,
		White,
		Gray;

		public int getValue()
		{
			return this.ordinal();
		}

		public static ColorGroup forValue(int value)
		{
			return values()[value];
		}
	}

	public static class Sku
	{

		public Sku(String code, String description, double large, double width, double height, double volume, double weight, String section, String type, String subsection, double cost, double price, int colorIndex, double ltFcst, double rollFcst)
		{
			this.setcode(code);
			this.setdescription(description);
			this.setlarge(large);
			this.setwidth(width);
			this.setheight(height);
			this.setvolume(volume);
			this.setweight(weight);
			this.settype(type);
			this.setsection(section);
			this.setsubsection(subsection);
			this.setcost(cost);
			this.setprice(price);
			SkuClustering skuClust = new SkuClustering();
			this.setcolor(skuClust.GetColor(colorIndex));
			this.setR(getcolor().R);
			this.setG(getcolor().G);
			this.setB(getcolor().B);
			this.setltFcst(ltFcst);
			this.setrollFcst(rollFcst);
		}

		private String code;
		public final String getcode()
		{
			return code;
		}
		public final void setcode(String value)
		{
			code = value;
		}
		private String description;
		public final String getdescription()
		{
			return description;
		}
		public final void setdescription(String value)
		{
			description = value;
		}
		private String type;
		public final String gettype()
		{
			return type;
		}
		public final void settype(String value)
		{
			type = value;
		}
		private double large;
		public final double getlarge()
		{
			return large;
		}
		public final void setlarge(double value)
		{
			large = value;
		}
		private double width;
		public final double getwidth()
		{
			return width;
		}
		public final void setwidth(double value)
		{
			width = value;
		}
		private double height;
		public final double getheight()
		{
			return height;
		}
		public final void setheight(double value)
		{
			height = value;
		}
		private double volume;
		public final double getvolume()
		{
			return volume;
		}
		public final void setvolume(double value)
		{
			volume = value;
		}
		private double weight;
		public final double getweight()
		{
			return weight;
		}
		public final void setweight(double value)
		{
			weight = value;
		}
		private String section;
		public final String getsection()
		{
			return section;
		}
		public final void setsection(String value)
		{
			section = value;
		}
		private String subsection;
		public final String getsubsection()
		{
			return subsection;
		}
		public final void setsubsection(String value)
		{
			subsection = value;
		}
		private double cost;
		public final double getcost()
		{
			return cost;
		}
		public final void setcost(double value)
		{
			cost = value;
		}
		private double price;
		public final double getprice()
		{
			return price;
		}
		public final void setprice(double value)
		{
			price = value;
		}
		private SkuClustering.Color color;
		public final SkuClustering.Color getcolor()
		{
			return color;
		}
		public final void setcolor(SkuClustering.Color value)
		{
			color = value;
		}
		private double R;
		public final double getR()
		{
			return R;
		}
		public final void setR(double value)
		{
			R = value;
		}
		private double G;
		public final double getG()
		{
			return G;
		}
		public final void setG(double value)
		{
			G = value;
		}
		private double B;
		public final double getB()
		{
			return B;
		}
		public final void setB(double value)
		{
			B = value;
		}
		private double ltFcst;
		public final double getltFcst()
		{
			return ltFcst;
		}
		public final void setltFcst(double value)
		{
			ltFcst = value;
		}
		private double rollFcst;
		public final double getrollFcst()
		{
			return rollFcst;
		}
		public final void setrollFcst(double value)
		{
			rollFcst = value;
		}
		private double similWeight;
		public final double getsimilWeight()
		{
			return similWeight;
		}
		public final void setsimilWeight(double value)
		{
			similWeight = value;
		}


		@Override
		public String toString()
		{
			return getcode() + "\t" + getdescription() + "\t" + gettype() + "\t" + getlarge() + "\t" + getwidth() + "\t" + getheight() + "\t" + getvolume() + "\t" + getweight() + "\t" + getsection() + "\t" + getsubsection() + "\t" + getcost() + "\t" + getprice() + "\t" + getcolor().name;
		}
	}

	//endregion
}