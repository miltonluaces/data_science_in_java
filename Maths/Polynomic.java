package NumCalc;

//region Imports

import ClassicalStat.*;
import java.util.*;

//endregion



/**  Class for polynomic operations 
*/
public class Polynomic
{

	//region Constructor

	/**  Constructor 
	*/
	public Polynomic()
	{
	}

	//endregion

	//region Polynomic Regression

	/** 
	 Regresion cuadratica
	 a0 + a1x + a2x2   -  Sistema : AU = Y
	 
	 @param alX independent variable 
	 @param alY dependent variable 
	 @param grade grade of the polynom 
	 @return  list of polynomial coefficents 
	*/
	public final ArrayList<Double> Regression(ArrayList<Double> alX, ArrayList<Double> alY, int grade)
	{
		double[][] datosA = new double[alX.size()][];
		for (int i = 0;i < alX.size();i++)
		{
			datosA[i] = new double[grade+1];
			datosA[i][0] = 1;
			datosA[i][1] = (double)alX.get(i);
			for (int g = 2;g <= grade;g++)
			{
				datosA[i][g] = Math.pow(datosA[i][1], g);
			}
		}
		double[][] datosY = new double[alY.size()][];
		for (int i = 0;i < alY.size();i++)
		{
			datosY[i] = new double[1];
			datosY[i][0] = (double)alY.get(i);
		}

		Matrix A = new Matrix(datosA);
		Matrix Y = new Matrix(datosY);

		IMatrix X = A.Solve(Y);
		ArrayList<Double> solucion = new ArrayList<Double>();
		for (int i = 0;i < X.Rows;i++)
		{
			solucion.add(X[i][0]);
		}
		return solucion;
	}

	/**  Get Regression values 
	 @param X independent variable 
	 @param coeffs list of coefficents 
	 @return  list of regression values 
	*/
	public final ArrayList<Double> GetRegValues(ArrayList<Double> X, ArrayList<Double> coeffs)
	{
		ArrayList<Double> regValues = new ArrayList<Double>();
		double y;
		for (double x : X)
		{
			y = 0.0;
			for (int i = 0;i < coeffs.size();i++)
			{
				y += coeffs.get(i) * Math.pow(x,i);
			}
			regValues.add(y);
		}
		return regValues;
	}


	/**  Get Regression value 
	 @param x independent variable 
	 @param coeffs list of coefficents 
	 @return  regression value 
	*/
	public final double GetRegValue(double x, ArrayList<Double> coeffs)
	{
		double y = 0.0;
		for (int i = 0;i < coeffs.size();i++)
		{
			y += coeffs.get(i) * Math.pow(x, i);
		}
		return y;
	}

	//endregion

	//region M Regresion 

	/*
	public Sttatic List<double> GetResiduals(List<double> X, List<double> Y, List<double> coeffs) {
	    if(X.Count != Y.Count) { throw new Exception("X and Y must have the same size.");  }

	    List<double> residuals = new List<double>();
	    List<double> regValues = GetRegValues(X, coeffs);
	    for(int i=0;i<Y.Count;i++) {
	        residuals.Add(y[i] - regValues[i]);
	    }
	    return residuals;
	}

	public List<double> BmReg(List<double> xList, List<double> yList, int iter, double bend) {

	    //conversion de listas a matrices
	    Matrix x = new Matrix(xList.Count, 1);
	    for(int i=0;i<xList.Count;i++) { x[0, i] = xList[i]; }
	    Matrix y = new Matrix(yList.Count, 1);
	    for(int i=0;i<yList.Count;i++) { y[0, i] = yList[i]; }
	    
	    //valores por defecto
	    if(iter == -1) { iter = 20; }
	    if(bend == -1) { bend = 2 * Math.Sqrt(x.Columns + 1)/ x.Rows; }

	    //calculo de coeficientes y residuo inicial�
	    List<double> initCoeffs = Regression(xList, yList, 1);
	    List<double> residuals = GetResiduals(xList, yList, initCoeffs);
	    
	    
	    //parametros
	    Matrix x1 = CBind(x, 1);
	    double nu = Math.Sqrt(1 - Hat(x1));
	    int low = x.Columns + 1;
	    double eps = 0.0001;
	    double maxAbs = 0.0;
	    double newResid, wt;
	    List<double> newCoeffs = new List<double>();
	    List<double> ev = new List<double>();
	    List<double> rov;

	    //bucle principal
	    for(int i=0;i<iter;i++) {
	        
	        ev.Add(Math.Abs(resid));
	        ev.Sort();

	        //calculo nuevo ajuste
	        List<double> scale = 0.0;
	        List<double> vals = new List<double>();
	        for(int i=low;i<=yList.Count;i++) { vals.Add(ev[i]); }
	        vals.Sort();
	        double median = vals[vals.Count/2];
	        scale = median/InvNormal_acum(0.75);  //scale = Median(ev[c(low:length(y))])/qnorm(0.75);
	        rov = (resid/scale)/nu; 
	        double psi = (Math.Abs(rov) <= bend)?  rov : bend * Math.Sign(rov);
	        List<double> wt = new List<double>();
	        for(int i=0;i<residuals.Count;i++) { wt.Add(nu * psi / (residuals[i]/scale)); }
	        newCoeffs = Regression(xList, yList, wt);
	        newResid = GetResiduals(xList, yList, newCoeffs);
	        
	        //criterio de parada: maxima diferencia < eps
	        double abs;
	        for(int j=0;j<initCoeffs.Count;j++) {
	            abs = Math.Abs(newCoeffs[j] - initCoeffs[j]);
	            if(abs > maxAbs) { maxAbs = abs; }
	        }
	        if(maxAbs < eps) { break; }
	    
	        //siguiente paso iterativo
	        initCoeffs = newCoeffs;
	        resid = newResid;
	    }

	    //obtenci�n del nuevo residuo y newCoeffs
	    IMatrix yMenosX1 = y.Subtraction(x1);
	    IMatrix newCoeffsMat = ToMatrix(newCoeffs);
	    IMatrix residMat = yMenosX1.Multiply(newCoeffsMat);
	    resid = residMat.Toscalar();
	   
	    if (maxAbs >= eps) { throw new Exception("failed to converge in" + iter + "steps"); }
	    
	    //coefficents	
	    return newCoeffs;
	}

	private double LsFitResiduals(List<double> x, List<double> y, List<double> coeffs) {
	    
	    double residuals = 0.0;
	    return residuals;
	}

	private List<double> LsFitCoeffs(Matrix x, Matrix y, double wt) {
	    List<double> coeffs = new List<double>();
	    return coeffs;
	}

	private List<double> LsFitCoeffs(Matrix x, Matrix y) {
	    List<double> coeffs = new List<double>();
	    return coeffs;
	}

	private List<double> LsFitCoeffs(Matrix y) {
	    List<double> coeffs = new List<double>();
	    return coeffs;
	}

	private double LsFitResiduals(Matrix x, Matrix y, double wt) {
	    double resid = 0.0;
	    return resid;
	}

	private double Hat(Matrix x) {
	    double hat = -1;
	    x = CBind(x, 1);
	    n = x.Rows;
	    IQrDecomposition q = x.GetQrDecomposition();
	    Matrix qr = q.UpperTriangularFactor.Addition(q.OrthogonalFactor);
	    x = qr;
	    //apply(qr.qy(x, diag(1, nrow = n, ncol = x$rank))^2, 1, sum)
	    return hat;
	}

	private double Qy(Matrix qrqr, Matrix y) { 
	    //if (!is.qr(qr)) stop("argument is not a QR decomposition")
	    
	    int n = qrqr.Rows;                  //n<- as.integer(nrow(qr$qr))
	    int p = qrqr.Columns;               //p <- as.integer(ncol(qr$qr))
	    int df = qrRank;                     //df <- as.integer(qr$rank)
	    double ny = (double)y.Columns;      //ny <- as.integer(NCOL(y)) storage.mode(y) <- "double"
	    if(y.Rows != n) { throw new Exception("qr and y must have the same number of rows"); }
	    
	    Matrix qy = new Matrix(y.Rows, y.Columns);
	    for(int i=0;i<y.Rows;i++) {
	        for(int j=0;j<y.Columns;j++) {
	            qy[i,j] = n * ny;
	        }
	    }

	    //.Fortran("dqrqy", as.double(qr$qr), n, df, as.double(qr$qraux),y, ny, qy = qy, PACKAGE = "base")$qy
	}
	
	private Matrix CBind(Matrix x, double val) {
	    Matrix binded = new Matrix(x.Rows, x.Columns + 1);
	    for(int i=0;i<x.Rows;i++) {
	        binded[i, 0] = x[i, 0];
	        binded[i, 1] = val;
	    }
	    return binded;
	}

	private IMatrix ToMatrix(List<double> x) {
	    IMatrix xMat = new Matrix(x.Count, 1);
	    for(int i=0;i<x.Count;i++) { xMat[i, 0] = x[i]; }
	    return xMat;
	}

	public Sttatic double InvNormal_acum(double p) {
	    double q, r;

	    if(p < 0 || p > 1) {
	        return 0.0;
	    } 
	    else if(p == 0) {
	        return 999999999;
	    } 
	    else if(p == 1) {
	        return 999999999;
	    }
	    else if(p < LOW) {
	        //Rational approximation for lower region 
	        q = Math.Sqrt(-2*Math.Log(p));
	        return (((((c[0]*q+c[1])*q+c[2])*q+c[3])*q+c[4])*q+c[5]) /
				((((d[0]*q+d[1])*q+d[2])*q+d[3])*q+1);
	    } 
	    else if(p > HIGH) {
	        //Rational approximation for upper region 
	        q  = Math.Sqrt(-2*Math.Log(1-p));
	        return -(((((c[0]*q+c[1])*q+c[2])*q+c[3])*q+c[4])*q+c[5]) /
				((((d[0]*q+d[1])*q+d[2])*q+d[3])*q+1);
	    } 
	    else {
	        //Rational approximation for central region 
	        q = p - 0.5;
	        r = q*q;
	        return (((((a[0]*r+a[1])*r+a[2])*r+a[3])*r+a[4])*r+a[5])*q /
				(((((b[0]*r+b[1])*r+b[2])*r+b[3])*r+b[4])*r+1);
	    }
	}
	*/

	//endregion

	//region Weighted regression


	/**  Weighted Least Square Reg Known weights 
	 @param X independent variable 
	 @param Y dependent variable 
	 @param W list of weights 
	 @return  Result object that contains coefficents 
	*/
	public final Result WLSRegression(ArrayList<Double> X, ArrayList<Double> Y, ArrayList<Double> W)
	{
		if (X.size() != Y.size())
		{
			throw new RuntimeException(Strings.X_and_Y_must_have_the_same_size);
		}
		if (X.size() != W.size())
		{
			throw new RuntimeException(Strings.X_and_W_must_have_the_same_size);
		}

		int n = X.size();
		Result res = new Result(n, W);

		int i;
		double ss = 0.0;
		double sx = 0.0;
		double sy = 0.0;
		double mx, t;
		double st2 = 0.0;

		for (i = 0;i < n;i++)
		{
			ss += res.W.get(i);
			sx += X.get(i) * res.W.get(i);
			sy += Y.get(i) * res.W.get(i);
		}
		mx = sx / ss;

		res.c1 = 0.0;

		for (i = 0;i < n;i++)
		{
			t = X.get(i) - mx;
			st2 += Math.pow(t,2) * res.W.get(i);
			res.c1 += t * Y.get(i) * res.W.get(i);
		}

		res.c1 /= st2;
		res.c0 = (sy - sx * res.c1) / ss;

		res.sig0 = Math.sqrt((1.0 + sx * sx / (n * st2)) / ss);
		res.sig1 = Math.sqrt(1.0 / st2);

		res.chi2 = 0.0;
		for (i = 0;i < n;i++)
		{
			t = Y.get(i) - res.c0 - res.c1 * X.get(i);
			res.chi2 += Math.pow(t,2) * res.W.get(i);
		}
		return res;
	}


	/**   Weighted Least Square Reg Unknown weights 
	 @param X independent variable time series 
	 @param Y dependent variable time series 
	 @return  Result of calculation with coefficents 
	*/
	public final Result WLSRegression(ArrayList<Double> X, ArrayList<Double> Y)
	{
		if (X.size() != Y.size())
		{
			throw new RuntimeException(Strings.X_and_Y_must_have_the_same_size);
		}

		int n = X.size();
		Result res = new Result(n);
		res = WLSRegression(X, Y, res.W);

		res.ResetW();

		double fit, result, err;
		double max = 0.0;

		//find the datum with the largest abs(residual) 
		for (int i = 0;i < n;i++)
		{
			fit = X.get(i) * res.c1 + res.c0;
			result = Y.get(i) - fit;
			err = Math.abs(result);
			res.W.set(i, err);
			if (res.W.get(i).compareTo(max) > 0)
			{
				max = res.W.get(i);
			}
		}

		//scale the weights so that the most distant outlier is ignored altogether, and other weights are proportional to error 
		for (int i = 0;i < n;i++)
		{
			res.W.set(i, 1.0 - res.W.get(i) / max);
		}

		//use the calculated weights
		return WLSRegression(X, Y, res.W);
	}

	/**  Weighted Least Square Reg Unknown weights only Y values 
	 @param Y list of values 
	 @return  result of calculation 
	*/
	public final Result WLSRegression(ArrayList<Double> Y)
	{
		ArrayList<Double> X = new ArrayList<Double>();
		for (int i = 0;i < Y.size();i++)
		{
			X.add((double)i);
		}
		return WLSRegression(X, Y);
	}


	/**  Least Squared Regression with unknown weights (only Y values) 
	 @param Y values for regression 
	 @return  result of calculation 
	*/
	public final Result LSRegression(ArrayList<Double> Y)
	{
		ArrayList<Double> X = new ArrayList<Double>();
		for (int i = 0;i < Y.size();i++)
		{
			X.add((double)i);
		}
		ArrayList<Double> coeffs = Regression(X, Y, 1);
		Result res = new Result(Y.size());
		res.c0 = coeffs.get(0);
		res.c1 = coeffs.get(1);
		return res;
	}

	//endregion

	//region Polynoms derivatives

	/**  Derive from a list of coefficents 
	 @param coeffs list of coefficents 
	 @return  derived polynom 
	*/
	private ArrayList<Double> Derive(ArrayList<Double> coeffs)
	{
		ArrayList<Double> derived = new ArrayList<Double>();
		for (int i = 1;i < coeffs.size();i++)
		{
			derived.add(i * coeffs.get(i));
		}
		return derived;
	}

	/**  Derive on any order from a list of coefficents 
	 @param coeffs list of coefficents 
	 @param order derivative order 
	 @return  derived polynom 
	*/
	public final ArrayList<Double> Derive(ArrayList<Double> coeffs, int order)
	{
		ArrayList<Double> derived = new ArrayList<Double>();
		for (int i = 0;i < order;i++)
		{
			coeffs = Derive(coeffs);
		}
		return coeffs;
	}

	/**  Derive from a list of values 
	 @param serie list of values 
	 @param order derivative order 
	 @return  derived polynom 
	*/
	public final ArrayList<Double> DeriveSerie(ArrayList<Double> serie, int order)
	{
		if (serie.size() == 1)
		{
			ArrayList<Double> res = new ArrayList<Double>();
			for (int i = 0;i < serie.size();i++)
			{
				res.add(0.0);
			}
			return res;
		}
		int grade = serie.size() - 1;
		if (grade >= 7)
		{
			grade = 6;
		}
		ArrayList<Double> x = new ArrayList<Double>();
		for (int i = 0;i < serie.size();i++)
		{
			x.add((double)i);
		}
		ArrayList<Double> coeffs = Regression(x, serie, grade);
		coeffs = Derive(coeffs, order);
		return GetRegValues(x, coeffs);
	}

	/**  Calculate indexes of points where a change from negative to positive takes place 
	 @param serie the time series 
	 @return  list of indexes 
	*/
	public final ArrayList<Integer> ZerosNegPosIndexes(ArrayList<Double> serie)
	{
		boolean neg = false;
		ArrayList<Integer> zerosIndexes = new ArrayList<Integer>();
		for (int i = 0;i < serie.size();i++)
		{
			if (serie.get(i).compareTo(0) >= 0)
			{
				if (neg == true && i < serie.size() - 1)
				{
					zerosIndexes.add(i + 1);
				}
				neg = false;
			}
			else
			{
				neg = true;
			}
		}
		return zerosIndexes;
	}

	/**  Calculate indexes of points where a change from positive to negative takes place 
	 @param serie the time series 
	 @return  list of indexes 
	*/
	public final ArrayList<Integer> ZerosPosNegIndexes(ArrayList<Double> serie)
	{
		boolean pos = false;
		ArrayList<Integer> zerosIndexes = new ArrayList<Integer>();
		for (int i = 0;i < serie.size();i++)
		{
			if (serie.get(i).compareTo(0) >= 0)
			{
				if (pos == true && i < serie.size() - 1)
				{
					zerosIndexes.add(i + 1);
				}
				pos = false;
			}
			else
			{
				pos = true;
			}
		}
		return zerosIndexes;
	}

	//endregion


}