package DeepLearning;

import java.util.*;

	/**  Collection of Neurons 
	*/

public class NeuronaCollection extends CollectionBase
{

	public NeuronaCollection()
	{

	}

	public final void agregar(Neurona ne)
	{
			this.InnerList.add(ne);
	}

	public final void borrar(Neurona ne)
	{
			if (this.InnerList.contains(ne))
			{
				int indice = this.InnerList.indexOf(ne);
				this.RemoveAt(indice);
			}
	}

	public final void ordenar()
	{
			Collections.sort(this.InnerList);
	}

	public final void agregarArray(Neurona[] array)
	{
			for (Neurona ne : array)
			{
				this.agregar(ne);
			}
	}

	public final Object[] getArray()
	{
			return this.InnerList.toArray(new Object[0]);
	}
}