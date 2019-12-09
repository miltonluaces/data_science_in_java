package Clustering;

import NumCalc.*;
import ClassicalStat.*;
import java.util.*;

//region Dimension class

/**  Dimension struct 
*/
public class Dimension
{
	/**  dimension index 
	*/
	public int index;
	/**  number of divisions for this dimension 
	*/
	public int nDivisions;
	/**  if it shoud be logarithmed 
	*/
	public boolean logarithm;
	/**  if it should be normalized 
	*/
	public boolean normalize;
	/**  weight to be applied to this dimension 
	*/
	public double weight;
	/**  if it should be sorted 
	*/
	public boolean sort;
	/**  minimum value 
	*/
	public double min;
	/**  maximum value 
	*/
	public double max;
}