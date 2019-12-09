package Clustering;

public class ClustCalculator
{

	//region Fields

	private Dato dato;
	private KDTree kdt;


	//endregion

	//region Constructor

	/**  Constructor 
	 @param kdt KDTree to use 
	 @param dato value 
	*/
	public ClustCalculator(KDTree kdt, Dato dato)
	{
		this.kdt = kdt;
		this.dato = dato;
	}

	//endregion

	//region Public Methods

	/**  Calculating 
	*/
	public final void Calculate()
	{
		Dato centro = kdt.NearestNeighbourSearch(dato).dato;
		kdt.GetCluster(centro).Add(dato);
	}

	//endregion

}