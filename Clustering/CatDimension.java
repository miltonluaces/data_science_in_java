package Clustering;

import NumCalc.*;
import ClassicalStat.*;
import java.util.*;

public class CatDimension
{
	public int index;
	public double weight;
	public int nCats;
	public String[] catNames;
	public boolean[] catValues;

	public final void SetCategory(String catName)
	{
		for (int i = 0; i < nCats; i++)
		{
			if (catNames[i].equals(catName))
			{
				catValues[i] = true;
			}
			else
			{
				catValues[i] = false;
			}
		}
	}

	public final void SetCategory(int index)
	{
		catValues = new boolean[nCats];
		catValues[index] = true;
	}

	public final int GetIndex(String catName)
	{
		for (int i = 0; i < catNames.length; i++)
		{
			if (catNames[i].equals(catName))
			{
				return i;
			}
		}
		return -1;
	}

	public final int[] ToIntArray()
	{
		int[] intCats = new int[nCats];
		for (int i = 0;i < nCats;i++)
		{
			if (catValues[i] == true)
			{
				intCats[i] = 1;
			}
			else
			{
				intCats[i] = 0;
			}
		}
		return intCats;
	}
}
//endregion

