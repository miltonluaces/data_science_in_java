package Clustering;

import java.util.*;

//region Imports


//endregion


public class TestData
{

	//region Fields

	//endregion

	//region Constructor

	public TestData()
	{
	}

	//endregion

	//region Public Methods

	public final ArrayList<SkuClustering.Sku> GetData()
	{
		ArrayList<SkuClustering.Sku> skus = new ArrayList<SkuClustering.Sku>();
		java.io.InputStreamReader sr = new java.io.InputStreamReader("C:\\ProyectosVisualStudio\\NewClustering\\Algoritmos\\ClustData.csv");
		String line = sr.ReadLine();
		SkuClustering.Sku sku;
		while ((line = sr.ReadLine()) != null)
		{
			String[] tokens = line.split("[\\t]", -1);
			String code = tokens[0];
			String description = tokens[1];
			double large = Double.parseDouble(tokens[2]);
			double width = Double.parseDouble(tokens[3]);
			double height = Double.parseDouble(tokens[4]);
			double volume = Double.parseDouble(tokens[5]);
			double weight = Double.parseDouble(tokens[6]);
			String type = tokens[7];
			String section = tokens[8];
			String subsection = tokens[9];
			double cost = Double.parseDouble(tokens[10]);
			double price = Double.parseDouble(tokens[11]);
			int colorCode = Integer.parseInt(tokens[13]);
			double ltFcst = Double.parseDouble(tokens[14]);
			double rollFcst = Double.parseDouble(tokens[15]);
			sku = new SkuClustering.Sku(code, description, large, width, height, volume, weight, section, type, subsection, cost, price, colorCode, ltFcst, rollFcst);
			skus.add(sku);
		}
		sr.close();
		return skus;
	}
}
	//endregion

