package DeepLearning;

import java.util.*;

//region Imports


//endregion 

public class NNActivation
{

	//region Fields

	@FunctionalInterface
	public interface Activation
	{
		void invoke();
	}
	private ArrayList<Activation> acts;
	private boolean backpropagate;

	//endregion

	//region Constructors
	public NNActivation()
	{
		this(true);
	}

	public NNActivation(boolean backprop)
	{
		this.acts = new ArrayList<Activation>();
		this.backpropagate = backprop;
	}

	//endregion

	//region Properties

	public final ArrayList<Activation> getActs()
	{
		return acts;
	}
	public final void setActs(ArrayList<Activation> value)
	{
		acts = value;
	}
	public final boolean getBackpropagate()
	{
		return backpropagate;
	}
	public final void setBackpropagate(boolean value)
	{
		backpropagate = value;
	}

	//endregion

	//region Public Methods

	//region Activation

	public final Neuron Activate(ActFunction actFun, Neuron m)
	{
		Neuron res = new Neuron(m.Rows, m.Cols);
		int n = m.W.length;
		for (int i = 0; i < n; i++)
		{
			res.W[i] = actFun.Forward(m.W[i]);
		}
		if (backpropagate)
		{
			Activation act = () ->
			{
					for (int i = 0; i < n; i++)
					{
						m.Dw[i] += actFun.Backward(m.W[i]) * res.Dw[i];
					}
			};
			acts.add(act);
		}
		return res;
	}

	public final void Backward()
	{
		for (int i = acts.size() - 1; i >= 0; i--)
		{
			acts.get(i)();
		}
	}

	//endregion

	//region Matrix Operations
	public final Neuron Concat(Neuron m1, Neuron m2)
	{
		if (m1.Cols > 1 || m2.Cols > 1)
		{
			throw new RuntimeException("Expected column vectors");
		}

		Neuron res = new Neuron(m1.Rows + m2.Rows);
		int loc = 0;
		for (int i = 0; i < m1.W.length; i++)
		{
			res.W[loc] = m1.W[i];
			res.Dw[loc] = m1.Dw[i];
			res.Cache[loc] = m1.Cache[i];
			loc++;
		}
		for (int i = 0; i < m2.W.length; i++)
		{
			res.W[loc] = m2.W[i];
			res.Dw[loc] = m2.Dw[i];
			res.Cache[loc] = m2.Cache[i];
			loc++;
		}
		if (backpropagate)
		{
			Activation act = () ->
			{
					int index0 = 0;
					for (int i = 0; i < m1.W.length; i++)
					{
						m1.W[i] = res.W[index0];
						m1.Dw[i] = res.Dw[index0];
						m1.Cache[i] = res.Cache[index0];
						index0++;
					}
					for (int i = 0; i < m2.W.length; i++)
					{
						m2.W[i] = res.W[index0];
						m2.Dw[i] = res.Dw[index0];
						m2.Cache[i] = res.Cache[index0];
						index0++;
					}
			};

			acts.add(act);
		}
		return res;
	}

	public final Neuron Add(Neuron m1, Neuron m2)
	{
		if (m1.Rows != m2.Rows || m1.Cols != m2.Cols)
		{
			throw new RuntimeException("matrix dimension mismatch");
		}

		Neuron res = new Neuron(m1.Rows, m1.Cols);
		for (int i = 0; i < m1.W.length; i++)
		{
			res.W[i] = m1.W[i] + m2.W[i];
		}
		if (backpropagate)
		{
			Activation act = () ->
			{
					for (int i = 0; i < m1.W.length; i++)
					{
						m1.Dw[i] += res.Dw[i];
						m2.Dw[i] += res.Dw[i];
					}
			};
			acts.add(act);
		}
		return res;
	}

	public final Neuron Prod(Neuron m1, Neuron m2)
	{
		if (m1.Cols != m2.Rows)
		{
			throw new RuntimeException("matrix dimension mismatch");
		}

		int m1Rows = m1.Rows;
		int m1Cols = m1.Cols;
		int m2Cols = m2.Cols;
		Neuron res = new Neuron(m1Rows, m2Cols);
		int outcols = m2Cols;
		for (int i = 0; i < m1Rows; i++)
		{
			int m1Col = m1Cols * i;
			for (int j = 0; j < m2Cols; j++)
			{
				double dot = 0;
				for (int k = 0; k < m1Cols; k++)
				{
					dot += m1.W[m1Col + k] * m2.W[m2Cols * k + j];
				}
				res.W[outcols * i + j] = dot;
			}
		}
		if (backpropagate)
		{
			Activation act = () ->
			{
					for (int i = 0; i < m1.Rows; i++)
					{
						int outcol = outcols * i;
						for (int j = 0; j < m2.Cols; j++)
						{
							double b = res.Dw[outcol + j];
							for (int k = 0; k < m1.Cols; k++)
							{
								m1.Dw[m1Cols * i + k] += m2.W[m2Cols * k + j] * b;
								m2.Dw[m2Cols * k + j] += m1.W[m1Cols * i + k] * b;
							}
						}
					}

			};
			acts.add(act);
		}
		return res;
	}

	public final Neuron OneMinus(Neuron m)
	{
		Neuron ones = Neuron.Rep(m.Rows, m.Cols, 1.0);
		return Subtract(ones, m);
	}

	public final Neuron Subtract(Neuron m1, Neuron m2)
	{
		return Add(m1, Neg(m2));
	}

	public final Neuron Smul(Neuron m, double s)
	{
		Neuron m2 = Neuron.Rep(m.Rows, m.Cols, s);
		return Elmul(m, m2);
	}

	public final Neuron Neg(Neuron m)
	{
		Neuron negOnes = Neuron.Rep(m.Rows, m.Cols, -1.0);
		return Elmul(negOnes, m);
	}

	public final Neuron Elmul(Neuron m1, Neuron m2)
	{
		if (m1.Rows != m2.Rows || m1.Cols != m2.Cols)
		{
			throw new RuntimeException("matrix dimension mismatch");
		}

		Neuron res = new Neuron(m1.Rows, m1.Cols);
		for (int i = 0; i < m1.W.length; i++)
		{
			res.W[i] = m1.W[i] * m2.W[i];
		}
		if (backpropagate)
		{
			Activation act = () ->
			{
					for (int i = 0; i < m1.W.length; i++)
					{
						m1.Dw[i] += m2.W[i] * res.Dw[i];
						m2.Dw[i] += m1.W[i] * res.Dw[i];
					}
			};
			acts.add(act);
		}
		return res;
	}

	//endregion

	//endregion

}