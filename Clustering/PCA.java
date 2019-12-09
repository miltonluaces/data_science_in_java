package Clustering;

//region Imports


//endregion


 public class PCA
 {

	//region Fields

	 private double[][] data;
	 private int nPoints;
	 private int nDim;
	 private double[] vars;
	 private double[][] scores;

	 //endregion

	//region Constructor

	 public PCA()
	 {
	 }

	//endregion

	//region Properties

	public final double[] getVars()
	{
		return vars;
	}

	public final double[][] getScores()
	{
		return scores;
	}

	//endregion

	//region Public Methods

	public final void Calculate(double[][] data, int nPoints, int nDim)
	{
		 this.data = data;
		 this.nPoints = nPoints;
		 this.nDim = nDim;
		 int info = 0;
		 vars = new double[0];
		 scores = new double[0][0];
		 tangible.RefObject<Integer> tempRef_info = new tangible.RefObject<Integer>(info);
		 tangible.RefObject<double[]> tempRef_vars = new tangible.RefObject<double[]>(vars);
		 tangible.RefObject<double[][]> tempRef_scores = new tangible.RefObject<double[][]>(scores);
		 Calculate(data, nPoints, nDim, tempRef_info, tempRef_vars, tempRef_scores);
	 scores = tempRef_scores.argValue;
	 vars = tempRef_vars.argValue;
	 info = tempRef_info.argValue;
	}

	public final double GetLTValue(int dataIndex, int dimIndex)
	{
		double lt = 0;
		for (int i = 0; i < nDim ; i++)
		{
			lt += scores[dimIndex][i] * data[dataIndex][i];
		}
		return lt;
	}

	public final double[][] GetLinearTransfs(int nSelDims)
	{
		double[][] lts = new double[data.length][nSelDims];
		for (int i = 0; i < lts.length; i++)
		{
			for (int j = 0; j < lts[0].length; j++)
			{
				lts[i][j] = GetLTValue(i, j);
			}
		}
		return lts;
	}

	public final double GetLTValue(double[] datum, int dimIndex)
	{
		double lt = 0;
		for (int i = 0; i < nDim; i++)
		{
			lt += scores[dimIndex][i] * datum[i];
		}
		return lt;
	}

	public final double[] GetLinearTransf(double[] datum, int nSelDims)
	{
		double[] lt = new double[nSelDims];
		for (int i = 0; i < nSelDims; i++)
		{
			lt[i] = GetLTValue(datum, i);
		}
		return lt;
	}

	//endregion

	//region Private Methods

	//Input: data : dataset, array[nData,nDim], nData: dataset size, nDim:   number of dimensions
	//Output: result  :  return code: -4, if SVD subroutine haven't converged --- -1, if wrong parameters has been passed (NPoints<0, NVars<1) --- 1, if task is solved
	//vars    :   array[nDim] variance values
	//vectors :   array[nDim, nDim] matrix, whose columns store basis vectors.
	private void Calculate(double[][] x, int npoints, int nvars, tangible.RefObject<Integer> info, tangible.RefObject<double[]> s2, tangible.RefObject<double[][]> v)
	{
		double[][] a = new double[0][0];
		double[][] u = new double[0][0];
		double[][] vt = new double[0][0];
		double[] m = new double[0];
		double[] t = new double[0];
		int i = 0;
		int j = 0;
		double mean = 0;
		double variance = 0;
		double skewness = 0;
		double kurtosis = 0;
		int i_ = 0;

		info.argValue = 0;
		s2.argValue = new double[0];
		v.argValue = new double[0][0];


		// Check input data
		if (npoints < 0 || nvars < 1)
		{
			info.argValue = -1;
			return;
		}
		info.argValue = 1;
		// Special case: NPoints=0
		if (npoints == 0)
		{
			s2.argValue = new double[nvars - 1 + 1];
			v.argValue = new double[nvars - 1 + 1][nvars - 1 + 1];
			for (i = 0; i <= nvars - 1; i++)
			{
				s2.argValue[i] = 0;
			}
			for (i = 0; i <= nvars - 1; i++)
			{
				for (j = 0; j <= nvars - 1; j++)
				{
					if (i == j)
					{
						v.argValue[i][j] = 1;
					}
					else
					{
						v.argValue[i][j] = 0;
					}
				}
			}
			return;
		}

		// Calculate means
		m = new double[nvars - 1 + 1];
		t = new double[npoints - 1 + 1];
		for (j = 0; j <= nvars - 1; j++)
		{
			for (i_ = 0; i_ <= npoints - 1; i_++)
			{
				t[i_] = x[i_][j];
			}
			tangible.RefObject<Double> tempRef_mean = new tangible.RefObject<Double>(mean);
			tangible.RefObject<Double> tempRef_variance = new tangible.RefObject<Double>(variance);
			tangible.RefObject<Double> tempRef_skewness = new tangible.RefObject<Double>(skewness);
			tangible.RefObject<Double> tempRef_kurtosis = new tangible.RefObject<Double>(kurtosis);
			samplemoments(t, npoints, tempRef_mean, tempRef_variance, tempRef_skewness, tempRef_kurtosis);
		kurtosis = tempRef_kurtosis.argValue;
		skewness = tempRef_skewness.argValue;
		variance = tempRef_variance.argValue;
		mean = tempRef_mean.argValue;
			m[j] = mean;
		}

		// Center, apply SVD, prepare output
		a = new double[Math.max(npoints, nvars) - 1 + 1][nvars - 1 + 1];

		for (i = 0; i <= npoints - 1; i++)
		{
			for (i_ = 0; i_ <= nvars - 1; i_++)
			{
				a[i][i_] = x[i][i_];
			}
			for (i_ = 0; i_ <= nvars - 1; i_++)
			{
				a[i][i_] = a[i][i_] - m[i_];
			}
		}

		for (i = npoints; i <= nvars - 1; i++)
		{
			for (j = 0; j <= nvars - 1; j++)
			{
				a[i][j] = 0;
			}
		}
		tangible.RefObject<Double> tempRef_u = new tangible.RefObject<Double>(u);
		tangible.RefObject<Double> tempRef_vt = new tangible.RefObject<Double>(vt);
		if (!svd.rmatrixsvd(a, Math.max(npoints, nvars), nvars, 0, 1, 2, s2, tempRef_u, tempRef_vt))
		{
		vt = tempRef_vt.argValue;
		u = tempRef_u.argValue;
			info.argValue = -4;
			return;
		}
	else
	{
		vt = tempRef_vt.argValue;
		u = tempRef_u.argValue;
	}
		if (npoints != 1)
		{
			for (i = 0; i <= nvars - 1; i++)
			{
				s2.argValue[i] = Math.sqrt(s2.argValue[i]) / (npoints - 1);
			}
		}

		v.argValue = new double[nvars - 1 + 1][nvars - 1 + 1];
		blas.copyandtranspose(vt, 0, nvars - 1, 0, nvars - 1, v, 0, nvars - 1, 0, nvars - 1);
	}

	private void BuildBasis(double[][] data, int nData, int nDim, tangible.RefObject<Integer> result, tangible.RefObject<double[]> vars, tangible.RefObject<double[][]> vectors)
	{
		result.argValue = 0;
		vars.argValue = new double[0];
		vectors.argValue = new double[0][0];
		Calculate(data, nData, nDim, result, vars, vectors);
		return;
	}

	public final void samplemoments(double[] x, int n, tangible.RefObject<Double> mean, tangible.RefObject<Double> variance, tangible.RefObject<Double> skewness, tangible.RefObject<Double> kurtosis)
	{
		int i = 0;
		double v = 0;
		double v1 = 0;
		double v2 = 0;
		double stddev = 0;

		mean.argValue = 0;
		variance.argValue = 0;
		skewness.argValue = 0;
		kurtosis.argValue = 0;


		// Initialization, (special case 'nData=0)'
		mean.argValue = 0;
		variance.argValue = 0;
		skewness.argValue = 0;
		kurtosis.argValue = 0;
		stddev = 0;
		if (n <= 0)
		{
			return;
		}

		// Mean
		for (i = 0; i <= n - 1; i++)
		{
			mean.argValue = mean.argValue + x[i];
		}
		mean.argValue = mean.argValue / n;

		// Variance
		if (n != 1)
		{
			v1 = 0;
			for (i = 0; i <= n - 1; i++)
			{
				v1 = v1 + Math.pow(x[i] - mean.argValue, 2);
			}
			v2 = 0;
			for (i = 0; i <= n - 1; i++)
			{
				v2 = v2 + (x[i] - mean.argValue);
			}
			v2 = Math.pow(v2, 2) / n;
			variance.argValue = (v1 - v2) / (n - 1);
			if ((double)(variance.argValue) < (double)(0))
			{
				variance.argValue = 0;
			}
			stddev = Math.sqrt(variance.argValue);
		}

		// Skewness and kurtosis
		if ((double)(stddev) != (double)(0))
		{
			for (i = 0; i <= n - 1; i++)
			{
				v = (x[i] - mean.argValue) / stddev;
				v2 = Math.pow(v, 2);
				skewness.argValue = skewness.argValue + v2 * v;
				kurtosis.argValue = kurtosis.argValue + Math.pow(v2, 2);
			}
			skewness.argValue = skewness.argValue / n;
			kurtosis.argValue = kurtosis.argValue / n - 3;
		}
	}

	//endregion
 }