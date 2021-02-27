package ClassicalStat;

//region Class Datum

/**  Data for hypothesis tests 
*/
public class Datum implements java.lang.Comparable
{

/**  value 
*/
public double value;
/**  sample to which belongs the datum 
*/
public int sample;
/**  calculated rank of this value 
*/
public double rank;

/**  Constructor 
@param value value of the datum 
@param sample sample to which belongs the datum  
*/
public Datum(double value, int sample)
{
this.value = value;
this.sample = sample;
}

//region IComparable Members

public final int compareTo(Object obj)
{
if (this.value < ((Datum)obj).value)
{
	return -1;
}
else if (this.value > ((Datum)obj).value)
{
	return 1;
}
else
{
	return 0;
}
}

//endregion
}