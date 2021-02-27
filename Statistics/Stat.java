package ClassicalStat;

import java.util.ArrayList;

//region Imports


//endregion


//DotNetHelpers_doc_comment_start  Class which provides classical statisticals 
//DotNetHelpers_doc_comment_end
public class Stat
{

    //region Required constants

    private final static double ITMAX = 100;
    private final static double EPS = (double) 3.0e-7;
    private final static double FPMIN = (double) 1.0e-30;
    //readonly static double LOW = (double) 02425;
    //readonly static double HIGH = (double ) 97575;
    private final static double M_LN_SQRT_2PI = (double) 0.918938533204672741780329736406;
    private final static double M_2PI = (double) 6.283185307179586476925286766559;
    // readonly static double C7 = (double) 4.67;
    // readonly static double C8 = (double) 6.66;
    // readonly static double C9 = (double) 6.73;
    // readonly static double C10 = (double) 13.32;

    private final static double S0 = 0.083333333333333333333; //  1/12
    private final static double S1 = 0.00277777777777777777778; //  1/360
    private final static double S2 = 0.00079365079365079365079365; //  1/1260
    private final static double S3 = 0.000595238095238095238095238; //  1/1680
    private final static double S4 = 0.0008417508417508417508417508; //  1/1188

    private final static double [] a = { -3.969683028665376e+01, 2.209460984245205e+02, -2.759285104469687e+02, 1.383577518672690e+02, -3.066479806614716e+01, 2.506628277459239e+00 };
    private final static double [] b = { -5.447609879822406e+01, 1.615858368580409e+02, -1.556989798598866e+02, 6.680131188771972e+01, -1.328068155288572e+01 };

    private final static double [] c = { -7.784894002430293e-03, -3.223964580411365e-01, -2.400758277161838e+00, -2.549732539343734e+00, 4.374664141464968e+00, 2.938163982698783e+00 };

    private final static double [] d = { 7.784695709041462e-03, 3.224671290700398e-01, 2.445134137142996e+00, 3.754408661907416e+00 };

    //endregion 

    //region Required Functions

    //DotNetHelpers_doc_comment_start 
    //DotNetHelpers_doc_comment_body   DESCRIPTION
    //DotNetHelpers_doc_comment_body     Evaluates the "deviance part"
    //DotNetHelpers_doc_comment_body 	 bd0(x,M) :=  M * D0(sample/M) = M*[ sample/M * log(sample/M) + 1 - (sample/M) ] =
    //DotNetHelpers_doc_comment_body 	           =  x * log(sample/M) + M - sample
    //DotNetHelpers_doc_comment_body     where M = E[X] = n*p (or = lambda), for	  x, M > 0
    //DotNetHelpers_doc_comment_body 
    //DotNetHelpers_doc_comment_body 	in a manner that should be stable (with small relative error)
    //DotNetHelpers_doc_comment_body 	for all x and np. In particular for sample/np close to 1, direct
    //DotNetHelpers_doc_comment_body     evaluation fails, and evaluation is based on the Taylor series
    //DotNetHelpers_doc_comment_body 																												*	of log((1+v)/(1-v)) with v = (sample-np)/(sample+np).
    //DotNetHelpers_doc_comment_body @param x Parametro x de la funci�n Deviance Part.
    //DotNetHelpers_doc_comment_body @param np Parametro np de la funci�n Deviance Part.
    //DotNetHelpers_doc_comment_body @return  result of calculation 
    //DotNetHelpers_doc_comment_end
    public static double bd0(double x, double np)
    {
        double ej, s, s1, v;
        int j;

        if (Math.abs(x-np) < 0.1*(x+np))
        {
            v = (x-np)/(x+np);
            s = (x-np)*v; // s using v -- change by MM
            ej = 2*x*v;
            v = v*v;
            for (j=1; ; j++)
            { // Taylor series
                ej *= v;
                s1 = s+ej/((j<<1)+1);
                if (s1==s) // last term was effectively 0
                {
                    return (s1);
                }
                s = s1;
            }
        }
        //DotNetHelpers multiline preserve/* else:  | x - np |  is not too small */
        return (x*Math.log(x/np)+np-x);
    }

    //DotNetHelpers_doc_comment_start 
    //DotNetHelpers_doc_comment_body   DESCRIPTION
    //DotNetHelpers_doc_comment_body
    //DotNetHelpers_doc_comment_body	   Computes the log of the error term in Stirling's formula.
    //DotNetHelpers_doc_comment_body      For n Greater Than 15, uses the series 1/12n - 1/360n^3 + ...
    //DotNetHelpers_doc_comment_body      For n Less or Equal than  15, integers or half-integers, uses stored values.
    //DotNetHelpers_doc_comment_body	     For other n Less than 15, uses lgamma directly (don't use this to
    //DotNetHelpers_doc_comment_body        write lgamma!)
    //DotNetHelpers_doc_comment_body
    //DotNetHelpers_doc_comment_body 																												*	of log((1+v)/(1-v)) with v = (sample-np)/(sample+np).
    //DotNetHelpers_doc_comment_body @param n Parametro n de la funci�n que devuelve el log del error de la formula de Stirling.
    //DotNetHelpers_doc_comment_body @return  result of calculation 
    //DotNetHelpers_doc_comment_end
    public static double stirlerr(double n)
    {
        //DotNetHelpers multiline preserve/*
        //DotNetHelpers multiline preserve  error for 0, 0.5, 1.0, 1.5, ..., 14.5, 15.0.
        //DotNetHelpers multiline preserve*/

        double [] sferr_halves = { 0.0, 0.1534264097200273452913848, 0.0810614667953272582196702, 0.0548141210519176538961390, 0.0413406959554092940938221, 0.03316287351993628748511048, 0.02767792568499833914878929, 0.02374616365629749597132920, 0.02079067210376509311152277, 0.01848845053267318523077934, 0.01664469118982119216319487, 0.01513497322191737887351255, 0.01387612882307074799874573, 0.01281046524292022692424986, 0.01189670994589177009505572, 0.01110455975820691732662991, 0.010411265261972096497478567, 0.009799416126158803298389475, 0.009255462182712732917728637, 0.008768700134139385462952823, 0.008330563433362871256469318, 0.007934114564314020547248100, 0.007573675487951840794972024, 0.007244554301320383179543912, 0.006942840107209529865664152, 0.006665247032707682442354394, 0.006408994188004207068439631, 0.006171712263039457647532867, 0.005951370112758847735624416, 0.005746216513010115682023589, 0.005554733551962801371038690 };
        double nn;

        if (n <= 15.0)
        {
            nn = n + n;
            if (nn == (int)nn)
            {
                return (sferr_halves[(int)nn]);
            }
            return (Stat.lngamma(n + 1.0) - (n + 0.5)* Math.log(n) + n - M_LN_SQRT_2PI);
        }

        nn = n*n;
        if (n>500)
        {
            return ((S0-S1/nn)/n);
        }
        if (n> 80)
        {
            return ((S0-(S1-S2/nn)/nn)/n);
        }
        if (n> 35)
        {
            return ((S0-(S1-(S2-S3/nn)/nn)/nn)/n);
        }
        //DotNetHelpers multiline preserve/* 15 < n <= 35 : */
        return ((S0-(S1-(S2-(S3-S4/nn)/nn)/nn)/nn)/n);
    }
    //endregion Required Functions

    //region Constructor

    //endregion

    //region Clasical statistical distribution function

    //DotNetHelpers_doc_comment_start 
    //DotNetHelpers_doc_comment_body Devuelve la funci�n Gamma (df).
    //DotNetHelpers_doc_comment_body Verificada con el sofware del proyecto r
    //DotNetHelpers_doc_comment_body La aproximaci�n se obtiene de : http://www.rskey.org/gamma.htm
    //DotNetHelpers_doc_comment_body 
    //DotNetHelpers_doc_comment_body @param k Valor a buscar para la funci�n Gamma.
    //DotNetHelpers_doc_comment_body @return  calculated value 
    //DotNetHelpers_doc_comment_end
    public static double gamma(double k)
    {
        double [] cof = {75122.6331530, 80916.6278952, 36308.2951477, 8687.24529705, 1168.92649479, 83.8676043424, 2.50662827511};
        double sum = (double) 0, prod = (double) 1, aux = 5.5 + k;

        for (int j = 0; j < 7; j++)
        {
            sum += cof[j] * Math.pow(k, j);
            prod *= (k + j);
        }
        return (sum / prod) * Math.pow(aux, k+0.5) * Math.pow(Math.E, -aux);
    }

    //
    // Devuelve la funci�n Ln (Gamma (df)))
    // Verificada con las tablas del Abramovitz y ademas con el software del proyecto R
    //

    //DotNetHelpers_doc_comment_start Devuelve el logar�tmo neperiano de la funci�n Gamma (df).
    //DotNetHelpers_doc_comment_body @param k Valor a buscar para la funci�n Gamma.
    //DotNetHelpers_doc_comment_body @return  calculated value 
    //DotNetHelpers_doc_comment_end
    public static double lngamma(double k)
    {
        double x = k, y = k, tmp, ser;
        double [] cof = {76.18009172947146, -86.50532032941677, 24.01409824083091, -1.2317395724500155, 0.1208650973866179e-2, -0.5395239384953e-5};

        tmp = x+5.5;
        tmp -= (x + 0.5) * Math.log(tmp);
        ser = 1.000000000190015;
        for (int j = 0; j <= 5; j++)
        {
            ser += cof[j]/++y;
        }
        return -tmp + Math.log(2.5066282746310005 * ser / x);
    }

    //DotNetHelpers_doc_comment_start Returns the beta B(z, w) function. Tested using tables from internet.
    //DotNetHelpers_doc_comment_body @param z z value of the beta function.
    //DotNetHelpers_doc_comment_body @param w w value of the beta function.
    //DotNetHelpers_doc_comment_body @return  calculated value 
    //DotNetHelpers_doc_comment_end
    public static double beta(double z, double w)
    {
        return Math.exp(lngamma(z) + lngamma(w) - lngamma(z + w));
    }

    //DotNetHelpers_doc_comment_start Returns the incomplete gamma P(a, x) function.  Tested using tables from internet.
    //DotNetHelpers_doc_comment_body @param a a value of the incomplete gamma P(a, x) function.
    //DotNetHelpers_doc_comment_body @param x x value of the incomplete gamma P(a, x) function.
    //DotNetHelpers_doc_comment_body @return  calculated value 
    //DotNetHelpers_doc_comment_end
    public static double gammap(double a, double x)
    {
        double gamser = 0.0, gammcf = 0.0, gln = 0.0;

        if ((x < 0.0) || (a <= 0.0))
        {
            throw new IllegalArgumentException("Strings.Invalid_arguments_the_incomplete_gamma_function_P_a__x");
        }

        if (x < (a + 1.0))
        {
            DotNetHelpers.RefObject<Double> tempRef_gamser = new DotNetHelpers.RefObject<Double>(gamser);
            DotNetHelpers.RefObject<Double> tempRef_gln = new DotNetHelpers.RefObject<Double>(gln);
            gser(tempRef_gamser, a, x, tempRef_gln);
        gln = tempRef_gln.argValue;
        gamser = tempRef_gamser.argValue;
            return gamser;
        }
        else
        {
            DotNetHelpers.RefObject<Double> tempRef_gammcf = new DotNetHelpers.RefObject<Double>(gammcf);
            DotNetHelpers.RefObject<Double> tempRef_gln2 = new DotNetHelpers.RefObject<Double>(gln);
            gcf(tempRef_gammcf, a, x, tempRef_gln2);
        gln = tempRef_gln2.argValue;
        gammcf = tempRef_gammcf.argValue;
            return 1.0 - gammcf;
        }
    }

    //DotNetHelpers_doc_comment_start Returns the incomplete gamma Q(a, x) = 1 - P(a, x) function.  Tested using tables from internet.
    //DotNetHelpers_doc_comment_body @param a a value of the incomplete gamma Q(a, x) function.
    //DotNetHelpers_doc_comment_body @param x x value of the incomplete gamma Q(a, x) function.
    //DotNetHelpers_doc_comment_body @return  calculated value 
    //DotNetHelpers_doc_comment_end
    public static double gammaq(double a, double x)
    {
        double gamser = 0.0, gammcf = 0.0, gln = 0.0;

        if ((x < 0.0) || (a <= 0.0))
        {
            throw new IllegalArgumentException("Strings.Invalid_arguments_the_incomplete_gamma_function_P_a__x");
        }

        if (x < (a + 1.0))
        {
            DotNetHelpers.RefObject<Double> tempRef_gamser = new DotNetHelpers.RefObject<Double>(gamser);
            DotNetHelpers.RefObject<Double> tempRef_gln = new DotNetHelpers.RefObject<Double>(gln);
            gser(tempRef_gamser, a, x, tempRef_gln);
        gln = tempRef_gln.argValue;
        gamser = tempRef_gamser.argValue;
            return 1.0 - gamser;
        }
        else
        {
            DotNetHelpers.RefObject<Double> tempRef_gammcf = new DotNetHelpers.RefObject<Double>(gammcf);
            DotNetHelpers.RefObject<Double> tempRef_gln2 = new DotNetHelpers.RefObject<Double>(gln);
            gcf(tempRef_gammcf, a, x, tempRef_gln2);
        gln = tempRef_gln2.argValue;
        gammcf = tempRef_gammcf.argValue;
            return gammcf;
        }
    }

    //DotNetHelpers_doc_comment_start Requiered function for gamma, beta and incomplete beta and gamma functions.
    //DotNetHelpers_doc_comment_body @param gamser gamser value of this function.
    //DotNetHelpers_doc_comment_body @param a a value of this function.
    //DotNetHelpers_doc_comment_body @param x x value of this function.
    //DotNetHelpers_doc_comment_body @param gln gln value of this function.
    //DotNetHelpers_doc_comment_end
    public static void gser(DotNetHelpers.RefObject<Double> gamser, double a, double x, DotNetHelpers.RefObject<Double> gln)
    {
        double sum, del, ap;

        gln.argValue = lngamma(a);

        if (x <= 0.0)
        {
            if (x < 0.0)
            {
                throw new IllegalArgumentException("Strings.Invalid_argument_for_function_gser");
            }
            gamser.argValue = 0.0;
            return;
        }
        else
        {
            ap = a;
            del = sum = 1.0 / a;
            for (int n = 1; n <= ITMAX; n++)
            {
                ++ap;
                del *= x / ap;
                sum += del;
                if (Math.abs(del) < Math.abs(sum) * EPS)
                {
                    gamser.argValue = sum * Math.exp(-x + a * Math.log(x) - gln.argValue);
                    return;
                }
            }
            throw new IllegalArgumentException("Strings.El_argumento_a_es_muy_grande_e_ITMAX_demasiado_peque�o_en_la_funcion_gser");
        }

    }

    //DotNetHelpers_doc_comment_start Requiered function for gamma, beta and incomplete beta and gamma functions.
    //DotNetHelpers_doc_comment_body @param gammcf gammcf value of this function.
    //DotNetHelpers_doc_comment_body @param a a value of this function.
    //DotNetHelpers_doc_comment_body @param x x value of this function.
    //DotNetHelpers_doc_comment_body @param gln gln value of this function.
    //DotNetHelpers_doc_comment_end
    public static void gcf(DotNetHelpers.RefObject<Double> gammcf, double a, double x, DotNetHelpers.RefObject<Double> gln)
    {
        double an, b, c, d, del, h;
        int i = 0;

        gln.argValue = lngamma(a);
        b = x + 1.0 - a;
        c = 1.0 / FPMIN;
        d = 1.0 / b;
        h = d;

        for (i = 1; i <= ITMAX; i++)
        {
            an = -i * (i - a);
            b += 2.0;
            d = an * d + b;
            if (Math.abs(d) < FPMIN)
            {
                d = FPMIN;
            }
            c = b + an / c;
            if (Math.abs(c) < FPMIN)
            {
                c = FPMIN;
            }
            d = 1.0 / d;
            del = c * d;
            h *= del;
            if (Math.abs(del - 1.0) < EPS)
            {
                break;
            }
        }

        if (i > ITMAX)
        {
            throw new IllegalArgumentException("Strings.El_argumento_a_es_muy_grande_e_ITMAX_demasiado_peque�o_en_la_funcion_gcf");
        }
        gammcf.argValue = Math.exp(-x + a * Math.log(x)- gln.argValue) * h;
    }

    //DotNetHelpers_doc_comment_start Returns the error function.  Tested using tables from internet.
    //DotNetHelpers_doc_comment_body @param x x value of the error.
    //DotNetHelpers_doc_comment_body @return  calculated value 
    //DotNetHelpers_doc_comment_end
    public static double erff(double x)
    {
        return x < 0.0 ? -gammap(0.5, Math.pow(x, 2)) : gammap(0.5, Math.pow(x, 2));
    }

    //DotNetHelpers_doc_comment_start Returns the complementary error function. Tested using tables from internet.
    //DotNetHelpers_doc_comment_body @param x x value of the error.
    //DotNetHelpers_doc_comment_body @return  calculated value 
    //DotNetHelpers_doc_comment_end
    public static double erffc(double x)
    {
        return x < 0.0 ? 1.0 + gammap(0.5, Math.pow(x, 2)) : gammaq(0.5, Math.pow(x, 2));
    }


    //DotNetHelpers_doc_comment_start 
    //DotNetHelpers_doc_comment_body Returns the probability distribution of Chi-Cuadrado P(Chi-Cuadrado | v)
    //DotNetHelpers_doc_comment_body Tested using internet page
    //DotNetHelpers_doc_comment_body http://members.aol.com/johnp71/pdfs.html
    //DotNetHelpers_doc_comment_body 
    //DotNetHelpers_doc_comment_body @param chi Value to be tested.
    //DotNetHelpers_doc_comment_body @param v Degrees of freedom of the distribution.
    //DotNetHelpers_doc_comment_body @return  calculated value 
    //DotNetHelpers_doc_comment_end
    public static double chi_square_p(double chi, double v)
    {
        return gammap(v / 2.0, chi / 2.0);
    }

    //DotNetHelpers_doc_comment_start 
    //DotNetHelpers_doc_comment_body Returns the complementary probability distribution of Chi-Cuadrado P(Chi-Cuadrado | v)
    //DotNetHelpers_doc_comment_body Tested using internet page
    //DotNetHelpers_doc_comment_body http://members.aol.com/johnp71/pdfs.html
    //DotNetHelpers_doc_comment_body 
    //DotNetHelpers_doc_comment_body @param chi Value to be tested.
    //DotNetHelpers_doc_comment_body @param v Degrees of freedom of the distribution.
    //DotNetHelpers_doc_comment_body @return  calculated value 
    //DotNetHelpers_doc_comment_end
    public static double chi_square_q(double chi, double v)
    {
        return gammaq(v / 2.0, chi / 2.0);
    }

    //DotNetHelpers_doc_comment_start Returns the incomplete Beta (a, b, x) = Ix(a, b) function. Tested using tables from internet.
    //DotNetHelpers_doc_comment_body @param a a Parameter of the incomplete Beta (a, b, x) function.
    //DotNetHelpers_doc_comment_body @param b b Parameter of the incomplete Beta (a, b, x) function.
    //DotNetHelpers_doc_comment_body @param x x Parameter of the incomplete Beta (a, b, x) function.
    //DotNetHelpers_doc_comment_body @return  calculated value 
    //DotNetHelpers_doc_comment_end
    public static double betai(double a, double b, double x)
    {
        double bt;

        if (x < 0.0 || x > 1.0)
        {
            throw new IllegalArgumentException("Strings.Invalid_argument_x_for_function_betai");
        }
        if (x == 0.0 || x == 1.0)
        {
            bt = 0.0;
        }
        else
        {
            bt = Math.exp(lngamma(a + b) - lngamma(a) -lngamma(b) + a * Math.log(x) + b * Math.log(1.0 - x));
        }
        if (x < (a + 1.0) / (a + b + 2.0))
        {
            return bt * betacf(a, b, x) / a;
        }
        else
        {
            return 1.0 - bt * betacf(b, a, 1.0 -x) / b;
        }
    }

    //DotNetHelpers_doc_comment_start Used by function betai. Evaluates continued fraction for incomplete beta function by modified Lentz`s method
    //DotNetHelpers_doc_comment_body @param a a Parameter of the incomplete Beta (a, b, x) function.
    //DotNetHelpers_doc_comment_body @param b b Parameter of the incomplete Beta (a, b, x) function.
    //DotNetHelpers_doc_comment_body @param x x Parameter of the incomplete Beta (a, b, x) function.
    //DotNetHelpers_doc_comment_body @return  calculated value 
    //DotNetHelpers_doc_comment_end
    public static double betacf(double a, double b, double x)
    {
        int m, m2;
        double aa, c, d, del, h, qab, qam, qap;

        qab = a + b;
        qap = a + 1.0;
        qam = a - 1.0;
        c = 1.0;
        d = 1.0 - qab * x / qap;
        if (Math.abs(d) < FPMIN)
        {
            d = FPMIN;
        }
        d = 1.0 / d;
        h = d;
        for (m = 1; m <= ITMAX; m++)
        {
            m2 = 2* m;
            aa = m * (b - m) * x / ((qam + m2) * (a + m2));
            d = 1.0 + aa * d;
            if (Math.abs(d) < FPMIN)
            {
                d= FPMIN;
            }
            c = 1.0 + aa / c;
            if (Math.abs(c) < FPMIN)
            {
                c= FPMIN;
            }
            d = 1.0 / d;
            h *= d * c;
            aa = -(a + m) * (qab + m) * x / ((aa + m2) * (qap + m2));
            d = 1.0 + aa * d;
            if (Math.abs(d) < FPMIN)
            {
                d= FPMIN;
            }
            c = 1.0 + aa / c;
            if (Math.abs(c) < FPMIN)
            {
                c= FPMIN;
            }
            d = 1.0 / d;
            del = d * c;
            h *= del;
            if (Math.abs(del - 1.0) < EPS)
            {
                break;
            }
        }
        if (m > ITMAX)
        {
            throw new IllegalArgumentException("Strings.Argument_a_or_b_are_too_big_or_MAXIT_is_too_small_in_function_betacf");
        }
        return h;
    }

    //
    // Devuelve la probabilidad x de un distribuci�n T-Student con n grados de libertad,
    // moda m y parametro de escala c, Tn(m, c)
    //

    //DotNetHelpers_doc_comment_start 
    //DotNetHelpers_doc_comment_body Returns the probability of a T-Student distribution with n degrees
    //DotNetHelpers_doc_comment_body of freedom, mode m and scale parameter c, Tn(m, c)
    //DotNetHelpers_doc_comment_body 
    //DotNetHelpers_doc_comment_body @param x x Value for which the probability is calculated .
    //DotNetHelpers_doc_comment_body @param m m Mode of the T-Student distribution.
    //DotNetHelpers_doc_comment_body @param c c Scale parameter of the T-Student distribution.
    //DotNetHelpers_doc_comment_body @param n n degrees of freedom of the T-Student distribution.
    //DotNetHelpers_doc_comment_body @param give_log give_log if true returns the log of the value, else just the value.
    //DotNetHelpers_doc_comment_body @return  calculated tStudent 
    //DotNetHelpers_doc_comment_end
    public static double TStudent(double x, double m, double c, double n, boolean give_log)
    {
//DotNetHelpers multiline preserve/*			return ((gamma ((n + 1) / 2) * Math.Pow (n, n / 2)) / 
//DotNetHelpers multiline preserve				(gamma (n / 2) * Math.Pow ((Math.PI * c), 0.5))) * 
//DotNetHelpers multiline preserve				Math.Pow ((n + Math.Pow (x - m, 2) / c), - (n + 1) / 2);
//DotNetHelpers multiline preserve*/
        double t, u, x_01 = (x - m) / c, x_01square = Math.pow(x_01, 2);

        t = -bd0(n/2.0, (n+1)/2.0) + stirlerr((n+1)/2.0) - stirlerr(n/2.0);
        if (x_01square > 0.2 * n)
        {
            u = Math.log(1 + x_01square / n) * n/2;
        }
        else
        {
            u = -bd0(n/2.0, (n+x_01square)/2.0) + Math.pow(x_01, 2)/2.0;
        }


        if (give_log)
        {
            return (-0.5*Math.log(M_2PI*(1+x_01square /n))+(t-u));
        }
        else
        {
            return (Math.exp(t-u) / Math.sqrt(M_2PI*(1+x_01square/n)));
        }
    }



    //DotNetHelpers_doc_comment_start 
    //DotNetHelpers_doc_comment_body Returns the log value of the probability of a T-Student distribution with n degrees
    //DotNetHelpers_doc_comment_body of freedom, mode m and scale parameter c, ln (Tn(m, c))
    //DotNetHelpers_doc_comment_body 
    //DotNetHelpers_doc_comment_body 0.5 * Ln (Pi * e) = 0.072364942924700087071713675676529
    //DotNetHelpers_doc_comment_body 
    //DotNetHelpers_doc_comment_body @param x x Value for which the probability is calculated .
    //DotNetHelpers_doc_comment_body @param m m Mode of the T-Student distribution.
    //DotNetHelpers_doc_comment_body @param c c Scale parameter of the T-Student distribution.
    //DotNetHelpers_doc_comment_body @param n n degrees of freedom of the T-Student distribution.
    //DotNetHelpers_doc_comment_body @return  calculated value 
    //DotNetHelpers_doc_comment_end
    public static double lnTStudent(double x, double m, double c, double n)
    {
        double n1 = n + 1;

        return Math.log(gamma(n1 / 2)) + n / 2 * Math.log(n) - Math.log(gamma(n / 2)) - 0.072364942924700087071713675676529 - n1 / 2 * Math.log(n + Math.pow(x - m, 2) / Math.E);
    }


    //DotNetHelpers_doc_comment_start 
    //DotNetHelpers_doc_comment_body Returns the inverse accumulated probability of a T-Student distribution with n degrees
    //DotNetHelpers_doc_comment_body of freedom, mode 0 and scale parameter 1, Tn(0, 1)
    //DotNetHelpers_doc_comment_body The number o degrees of freedom has to be greater than 1 (n>=1), 
    //DotNetHelpers_doc_comment_body if not the value is not correct.
    //DotNetHelpers_doc_comment_body https://svn.r-project.org/R/trunk/src/nmath/dt.c
    //DotNetHelpers_doc_comment_body 
    //DotNetHelpers_doc_comment_body @param p p Probability to search for
    //DotNetHelpers_doc_comment_body @param ndf ndf Degrees of freedom
    //DotNetHelpers_doc_comment_body @param lower_tail lower_tail True = gets (-inf, x). False = gets (sample, +inf).
    //DotNetHelpers_doc_comment_body @return  calculated quantil 
    //DotNetHelpers_doc_comment_end
    public static double TStudent_quantil(double p, double ndf, boolean lower_tail)
    {
        double eps=1e-12, P, prob, q, y, a, b, c, d, x;
        double M_PI_2=1.570796326794896619231321691640; // pi/2
        boolean neg;

        // Si el numero de grados de libertad es menor que 1, el valor es incorrecto
        // FIXME: This test should depend on  ndf  AND p  !!
        //        and in fact should be replaced by
        //        something like Abramowitz & Stegun 26.7.5 (p.949)
        // Se ha comprobado que la aproximaci�n de Abramowitz & Stegun 26.7.5 (p.949)
        // no funciona cuando el numero de grados de libertad es < 1
        // La aproximaci�n en cuesti�n es la que tenemos a continuaci�n

//DotNetHelpers multiline preserve/*			if(p<=0 || p>=1) return -1;
//DotNetHelpers multiline preserve
//DotNetHelpers multiline preserve			if (ndf < 1)
//DotNetHelpers multiline preserve			{
//DotNetHelpers multiline preserve				double x1,g1,g2,g3,g4;
//DotNetHelpers multiline preserve				x1=Stat.InvNormal_acum (p);
//DotNetHelpers multiline preserve				g1=(Math.Pow(x1,3.0)+x1)/4.0;
//DotNetHelpers multiline preserve				g2=(5.0*Math.Pow(x1,5.0)+16.0*Math.Pow(x1,3.0)+3.0*x1)/
//DotNetHelpers multiline preserve					96.0;
//DotNetHelpers multiline preserve				g3=(3.0*Math.Pow(x1,7.0)+19.0*Math.Pow(x1,5.0)+
//DotNetHelpers multiline preserve					17.0*Math.Pow(x1,3.0)-15.0*x1)/384.0;
//DotNetHelpers multiline preserve				g4=(79.0*Math.Pow(x1,9.0)+776.0*Math.Pow(x1,7.0)+
//DotNetHelpers multiline preserve					1482.0*Math.Pow(x1,5.0)-1920.0*Math.Pow(x1,3.0)-
//DotNetHelpers multiline preserve					945.0*x1)/92160.0;
//DotNetHelpers multiline preserve				return x1+g1/ndf+g2/Math.Pow(ndf,2.0)+
//DotNetHelpers multiline preserve					g3/Math.Pow(ndf,3.0)+g4/Math.Pow(ndf,4.0);
//DotNetHelpers multiline preserve			}
//DotNetHelpers multiline preserve*/
        if(p<=0 || p>=1 || ndf<1)
        {
            return -1;
        }
        if((lower_tail && p > 0.5) || (!lower_tail && p < 0.5))
        {
            neg = false;
            P = 2 * (lower_tail ? (1 - p) : p);
        }
        else
        {
            neg = true;
            P = 2 * (lower_tail ? p : (1 - p));
        }

        if(Math.abs(ndf - 2) < eps)
        { // df ~= 2
            q=Math.sqrt(2 / (P * (2 - P)) - 2);
        }
        else if (ndf < 1 + eps)
        { // df ~= 1
            prob = P * M_PI_2;
            q = Math.cos(prob)/Math.sin(prob);
        }
        else
        { //-- usual case;  including, e.g.,  df = 1.1
            a = 1 / (ndf - 0.5);
            b = 48 / (a * a);
            c = ((20700 * a / b - 98) * a - 16) * a + 96.36;
            d = ((94.5 / (b + c) - 3) / b + 1) * Math.sqrt(a * M_PI_2) * ndf;
            y = Math.pow(d * P, 2 / ndf);
            if (y > 0.05 + a)
            {
                //DotNetHelpers multiline preserve/* Asymptotic inverse expansion about normal */
                //x = qnorm(0.5 * P, false);
                x = R.qnorm(0.5 * P, 0, 1, false, false);
                y = x * x;
                if (ndf < 5)
                {
                    c += 0.3 * (ndf - 4.5) * (x + 0.6);
                }
                c = (((0.05 * d * x - 5) * x - 7) * x - 2) * x + b + c;
                y = (((((0.4 * y + 6.3) * y + 36) * y + 94.5) / c - y - 3) / b + 1) * x;
                y = a * y * y;
                if (y > 0.002)
                {
                    y = Math.exp(y) - 1;
                }
                else
                { // Taylor of    e^y -1 :
                    y = (0.5 * y + 1) * y;
                }
            }
            else
            {
                y = ((1 / (((ndf + 6) / (ndf * y) - 0.089 * d - 0.822) * (ndf + 2) * 3) + 0.5 / (ndf + 4)) * y - 1) * (ndf + 1) / (ndf + 2) + 1 / y;
            }
            q = Math.sqrt(ndf * y);
        }
        if(neg)
        {
            q = -q;
        }
        return q;
    }

    //DotNetHelpers_doc_comment_start 
    //DotNetHelpers_doc_comment_body Returns the distribution probability function T-Student A(t | v), with n degrees
    //DotNetHelpers_doc_comment_body of freedom, mode m and scale parameter c, Tn(m, c)
    //DotNetHelpers_doc_comment_body 
    //DotNetHelpers_doc_comment_body @param t t Value to be tested.
    //DotNetHelpers_doc_comment_body @param v v degrees of freedom of the T-Student distribution.
    //DotNetHelpers_doc_comment_body @param m m Mode of the T-Student distribution.
    //DotNetHelpers_doc_comment_body @param e e Scale parameter of the T-Student distribution.
    //DotNetHelpers_doc_comment_body @return  calculated tStrudent 
    //DotNetHelpers_doc_comment_end
    public static double Test_TStudent(double t, double v, double m, double e)
    {
        return 1.0 - betai(v / 2.0, 0.5, v / (v + Math.pow((t - m) / e, 2)));
    }

    //DotNetHelpers_doc_comment_start 
    //DotNetHelpers_doc_comment_body Returns the probability that c chi-2distrubution with dof degrees of freedom is
    //DotNetHelpers_doc_comment_body less or equal to x
    //DotNetHelpers_doc_comment_body 
    //DotNetHelpers_doc_comment_body @param x sample Value to be tested.
    //DotNetHelpers_doc_comment_body @param dof dof degrees of freedom of the T-Student distribution.
    //DotNetHelpers_doc_comment_body @return  calculated chi2 
    //DotNetHelpers_doc_comment_end
    public static double Test_Chi2(double x, double dof)
    {
        return 1 - (gammap(dof/2, x/2) / gamma(dof/2));
    }

    //DotNetHelpers_doc_comment_start 
    //DotNetHelpers_doc_comment_body Returns the probability function F, Q(F | v1, v2)
    //DotNetHelpers_doc_comment_body where t es the value to be tested and v1 and v2 the variances or degrees of freedom
    //DotNetHelpers_doc_comment_body Tested though url  http://members.aol.com/johnp71/pdfs.html
    //DotNetHelpers_doc_comment_body 
    //DotNetHelpers_doc_comment_body @param F F Value to be tested.
    //DotNetHelpers_doc_comment_body @param v1 v1 First variance or degrees of freedom of the F function.
    //DotNetHelpers_doc_comment_body @param v2 v2 Second variance or degrees of freedom of the F function.
    //DotNetHelpers_doc_comment_body @return  calculated value 
    //DotNetHelpers_doc_comment_end
    public static double Test_F(double F, double v1, double v2)
    {
        return 1 - betai(v2 / 2.0, v1 / 2.0, v2 / (v2 + (v1 * F)));
    }

    //DotNetHelpers_doc_comment_start  F Distribution 
    //DotNetHelpers_doc_comment_body @param f f statistic 
    //DotNetHelpers_doc_comment_body @param v1 degrees of freedom 1 
    //DotNetHelpers_doc_comment_body @param v2 degrees of freedom 2 
    //DotNetHelpers_doc_comment_body @return  f distribution value 
    //DotNetHelpers_doc_comment_end
    public static double FDistribution(double f, double v1, double v2)
    {
        if((v1*f)/(v1 + v2*f) > 1)
        {
            return 0;
        } //TODO: Verificar
        return betai(v1/2, v2/2, v1/(v1 + v2*f));
    }

    // the incomplete beta function from 0 to x with parameters a, b. x must be in (0,1) 
    private static double Betai(double a, double b, double x)
    {
        double bt=0;
        double beta = Double.MAX_VALUE;
        if(x==0 || x==1)
        {
            bt = 0;
        }
        else if((x>0)&&(x<1))
        {
            bt = gamma(a+b)*Math.pow(x, a)*Math.pow(1-x, b)/(gamma(a)*gamma(b));
        }
        if(x<(a+1)/(a+b+2))
        {
            beta = bt*betacf(a, b, x)/a;
        }
        else
        {
            beta = 1-bt*betacf(b, a, 1-x)/b;
        }
        return beta;
    }


    //endregion Clasical statistical distribution function

  //region Descriptive statistical function

  		/** 
  		 Returns the Mean value "Values".
  		 If FirstP is less than 0, returns -1
  		 If FirstP is greater than NValues, returns -2
  		 FirstP starts to count on 0.
  		 
  		 @param Values Values = Array of values from which to calculate the mean.
  		 @param NValues NValues = Number of values passed in the "Values" array.
  		 @param FirstP FirstP = First period to consider of the values, 0 all periods.
  		 @return  mean value 
  		*/
  		public static double Mean(double[] Values, int NValues, int FirstP)
  		{
  			double Sum = (double) 0;

  			if (FirstP < 0)
  			{
  				return (-1);
  			}
  			if (FirstP > NValues)
  			{
  				return (-2);
  			}

  			if (NValues != 0)
  			{
  				for (int i = FirstP; i < NValues; i++)
  				{
  					Sum += Values [i];
  				}
  				return (Sum / (NValues - FirstP));
  			}
  			else
  			{
  				return ((double) 0);
  			}
  		}

  		/** 
  		 Returns the Mean value "Values".
  		 
  		 @param Values Values = Array of values from which to calculate the mean.
  		 @return  mean values 
  		*/
  		public static double Mean(double[] Values)
  		{
  			double Sum = 0.0;
  			int NValues = Values.length;

  			if (NValues != 0)
  			{
  				for (int i = 0; i < NValues; i++)
  				{
  					Sum += Values[i];
  				}
  				return (Sum / NValues);
  			}
  			else
  			{
  				return (0.0);
  			}
  		}

  		/** 
  		 Returns the Mean value "Values".
  		 If FirstP is less than 0, returns -1
  		 If FirstP is greater than NValues, returns -2
  		 FirstP starts to count on 0.
  		 
  		 @param Values Values = Array of values from hich to calculate the mean.
  		 @param NValues NValues = Number of values passed in the "Values" array.
  		 @param UseValue UseValue = if 1, the value is used. If 0, the values is not used.
  		 @param FirstP FirstP = First period to consider of the values, 0 all periods.
  		 @return  calculated mean 
  		*/
  		public static double Mean(double[] Values, int NValues, int[] UseValue, int FirstP)
  		{
  			double Sum = (double) 0;
  			int values = 0;

  			if (FirstP < 0)
  			{
  				return (-1);
  			}
  			if (FirstP > NValues)
  			{
  				return (-2);
  			}

  			if (NValues != 0)
  			{
  				for (int i = FirstP; i < NValues; i++)
  				{
  					if (UseValue [i] == 1)
  					{
  						Sum += Values [i];
  						values++;
  					}
  				}
  				return (Sum / values);
  			}
  			else
  			{
  				return ((double) 0);
  			}
  		}

  		/** 
  		 Returns the varianza Muestral of array "Values".
  		 If FirstP is less than 0, returns -1
  		 If FirstP is greater than NValues, returns -2
  		 FirstP starts to count on 0.
  		 
  		 @param Values Values = Array of values from hich to calculate the mean.
  		 @param NValues NValues = Number of values passed in the "Values" array.
  		 @param FirstP FirstP = First period to consider of the values, 0 all periods.
  		 @return  variance value 
  		*/
  		public static double VarMuestral(double[] Values, int NValues, int FirstP)
  		{
  			if (Values.length == 1)
  			{
  				return 0;
  			}
  			int values = 0;
  			double Sum = (double) 0, SumSquares = (double) 0;

  			if (FirstP < 0)
  			{
  				return (-1);
  			}
  			if (FirstP > NValues)
  			{
  				return (-2);
  			}

  			for (int i = FirstP; i < NValues; i++)
  			{
  				Sum += Values [i];
  				SumSquares += Math.pow(Values [i], 2);
  				values++;
  			}
  			return (SumSquares - Math.pow(Sum, 2)) / (values - 1);
  		}

  		/** 
  		 Returns the varianza Muestral of array "Values".
  		 If FirstP is less than 0, returns -1
  		 If FirstP is greater than NValues, returns -2
  		 FirstP starts to count on 0.
  		 
  		 @param Values Values = Array of values from hich to calculate the mean.
  		 @return  variance value 
  		*/
  		public static double VarMuestral(double[] Values)
  		{
  			int NValues = Values.length;
  			if (NValues == 1)
  			{
  				return 0;
  			}
  			double Sum = 0.0, SumSquares = 0.0;

  			for (int i = 0; i < NValues; i++)
  			{
  				Sum += Values[i];
  				SumSquares += Math.pow(Values[i], 2);
  			}
  			return ((NValues * SumSquares - Math.pow(Sum, 2)) / (NValues * (NValues - 1)));
  		}

  		/** 
  		 Returns the varianza Muestral of array "Values".
  		 If FirstP is less than 0, returns -1
  		 If FirstP is greater than NValues, returns -2
  		 FirstP starts to count on 0.
  		 
  		 @param Values Values = Array of values from hich to calculate the mean.
  		 @param NValues NValues = Number of values passed in the "Values" array.
  		 @param UseValue UseValue = if 1, the value is used. If 0, the values is not used.
  		 @param FirstP FirstP = First period to consider of the values, 0 all periods.
  		 @return  calculated variance 
  		*/
  		public static double VarMuestral(double[] Values, int NValues, int[] UseValue, int FirstP)
  		{
  			int values = 0;
  			double Sum = (double) 0, SumSquares = (double) 0;

  			if (FirstP < 0)
  			{
  				return (-1);
  			}
  			if (FirstP > NValues)
  			{
  				return (-2);
  			}

  			for (int i = FirstP; i < NValues; i++)
  			{
  				if (UseValue [i] == 1)
  				{
  					Sum += Values [i];
  					SumSquares += Math.pow(Values [i], 2);
  					values++;
  				}
  			}
  			return ((values * SumSquares - Math.pow(Sum, 2)) / (values * (values - 1)));
  		}

  		/** 
  		 Returns the varianza Muestral of array "Values1" and "Values2".
  		 If FirstP is less than 0, returns -1
  		 If FirstP is greater than NValues, returns -2
  		 FirstP starts to count on 0.
  		 
  		 @param Values1 Values1 = First Array of values from which to calculate the mean.
  		 @param Values2 Values1 = Second Array of values from which to calculate the mean.
  		 @param NValues NValues = Number of values passed in the "Values" array.
  		 @param FirstP FirstP = First period to consider of the values, 0 all periods.
  		 @return  co variance value 
  		*/
  		public static double CoVarMuestral(double[] Values1, double[] Values2, int NValues, int FirstP)
  		{
  			double Mean1 = Mean(Values1, NValues, FirstP), Mean2 = Mean(Values2, NValues, FirstP);
  			double Sum = (double) 0;

  			if (FirstP < 0)
  			{
  				return (-1);
  			}
  			if (FirstP > NValues)
  			{
  				return (-2);
  			}

  			for (int i = FirstP; i < NValues; i++)
  			{
  				Sum += (Values1 [i] - Mean1) * (Values2 [i] - Mean2);
  			}
  			return (Sum / (NValues - FirstP));
  		}

  		/** 
  		 Returns the varianza Muestral of array "Values1" and "Values2".
  		 
  		 @param Values1 Values1 = First Array of values from which to calculate the mean.
  		 @param Values2 Values1 = Second Array of values from which to calculate the mean.
  		 @return  co variance value 
  		*/
  		public static double CoVarMuestral(double[] Values1, double[] Values2)
  		{
  			if (Values1.length != Values2.length)
  			{
  				throw new IllegalArgumentException("Strings.Both_vectors_must_have_the_same_number_of_values");
  			}

  			double Mean1 = Mean(Values1), Mean2 = Mean(Values2), Sum = (double)0;
  			int NValues = Values1.length;

  			for (int i = 0; i < NValues; i++)
  			{
  				Sum += (Values1[i] - Mean1) * (Values2[i] - Mean2);
  			}
  			return (Sum / NValues);
  		}

  		/** 
  		 Returns the varianza Muestral of array "Values1" and "Values2".
  		 If FirstP is less than 0, returns -1
  		 If FirstP is greater than NValues, returns -2
  		 FirstP starts to count on 0.
  		 
  		 @param Values1 Values1 = First Array of values from which to calculate the mean.
  		 @param Values2 Values2 = Second Array of values from which to calculate the mean.
  		 @param NValues NValues = Number of values passed in the "Values" array.
  		 @param UseValue UseValue = if 1, the value is used. If 0, the values is not used.
  		 @param FirstP FirstP = First period to consider of the values, 0 all periods.
  		 @return  co variance value 
  		*/
  		public static double CoVarMuestral(double[] Values1, double[] Values2, int NValues, int[] UseValue, int FirstP)
  		{
  			double Mean1 = Mean(Values1, NValues, UseValue, FirstP), Mean2 = Mean(Values2, NValues, UseValue, FirstP);
  			double Sum = (double) 0;
  			int values = 0;

  			if (FirstP < 0)
  			{
  				return (-1);
  			}
  			if (FirstP > NValues)
  			{
  				return (-2);
  			}

  			for (int i = FirstP; i < NValues; i++)
  			{
  				if (UseValue [i] == 1)
  				{
  					Sum += (Values1 [i] - Mean1) * (Values2 [i] - Mean2);
  					values++;
  				}
  			}
  			return (Sum / values);
  		}

  		/** 
  		 Returns the Correlation COeficient of array "Values1" and "Values2".
  		 If FirstP is less than 0, returns -1
  		 If FirstP is greater than NValues, returns -2
  		 FirstP starts to count on 0.
  		 
  		 @param Values1 Values1 = First Array of values from which to calculate the correlation coeficient.
  		 @param Values2 Values2 = Second Array of values from which to calculate the correlation coeficient.
  		 @param NValues NValues = Number of values passed in the "Values" array.
  		 @param UseValue UseValue = if 1, the value is used. If 0, the values is not used.
  		 @param FirstP FirstP = First period to consider of the values, 0 all periods.
  		 @return  correlation coefficent value 
  		*/
  		public static double CoefCorrelation(double[] Values1, double[] Values2, int NValues, int[] UseValue, int FirstP)
  		{
  			int NVal = 0;
  			double DNVal;

  			if (FirstP < 0)
  			{
  				return (-1);
  			}
  			if (FirstP > NValues)
  			{
  				return (-2);
  			}

  			for (int i = FirstP; i < NValues; i++)
  			{
  				if (UseValue [i] == 1)
  				{
  					NVal++;
  				}
  			}
  			DNVal = (double) NVal;
  			return ((DNVal / (DNVal - (double) 1.0)) * Stat.CoVarMuestral(Values1, Values2, NValues, UseValue, FirstP) / Math.sqrt(Stat.VarMuestral(Values1, NValues, UseValue, FirstP) * Stat.VarMuestral(Values2, NValues, UseValue, FirstP)));
  		}

  		/** 
  		 Returns the Correlation COeficient of array "Values1" and "Values2".
  		 
  		 @param Values1 Values1 = First Array of values from which to calculate the correlation coeficient.
  		 @param Values2 Values2 = Second Array of values from which to calculate the correlation coeficient.
  		 @return  correlation coefficent value 
  		*/
  		public static double CoefCorrelation(double[] Values1, double[] Values2)
  		{
  			if (Values1.length != Values2.length)
  			{
  				String msg = String.format("Strings.Both_Vectors_must_have_the_same_length, Values1.length, Values2.length");
  				throw new IllegalArgumentException(msg);
  			}

  			double DNVal = (double)Values1.length;
  			return ((DNVal / (DNVal - 1.0)) * Stat.CoVarMuestral(Values1, Values2) / Math.sqrt(Stat.VarMuestral(Values1) * Stat.VarMuestral(Values2)));
  		}

  		/** 
  		 Returns if Values1 is correlated with Values2, being Value the % of 
  		 probability to achieve.
  		 If FirstP is less than 0, Exception detected
  		 If FirstP is greater than NValues, Exception detected
  		 FirstP starts to count on 0.
  		 
  		 @param Values1 Values1c= First Array of values with which calculate correlation.
  		 @param Values2 Values2 = Second Array of values with which calculate correlation.
  		 @param NValues NValues = Number of values passed in the "Values" arrays.
  		 @param UseValue UseValue = if 1, the value is used. If 0, the values is not used.
  		 @param MinProb Percentage of probability with which the regressor must be selected.
  		 @param FirstP FirstP = First period to consider of the values, 0 all periods.
  		 @return  if the series are correlated 
  		*/
  		public static boolean Test_Correlation(double[] Values1, double[] Values2, int NValues, int[] UseValue, double MinProb, int FirstP)
  		{
  			double coefcorr, contraststat, limit;
  			int nperiods = 0, NObs;
  			if (FirstP < 0)
  			{
  				throw new IllegalArgumentException("Strings.Starting_Period_cannot_be_less_than_0");
  			}

  			if (FirstP > NValues)
  			{
  				throw new IllegalArgumentException("Strings.Starting_Period_cannot_greater_than_the_number_of_periods_in_the_Time_Series");
  			}

  			NObs = NValues - FirstP;
  			for (int i = FirstP; i < NValues; i++)
  			{
  				if (UseValue [i] == 1)
  				{
  					nperiods++;
  				}
  			}

  			limit = Stat.TStudent_quantil(MinProb / 100.0, nperiods - 2, true);
  			coefcorr = CoefCorrelation(Values1, Values2, NValues, UseValue, FirstP);
  			contraststat = coefcorr / Math.sqrt((1.0 - Math.pow(coefcorr, 2.0)) / (NValues - 2.0));
  			if (Double.isNaN(contraststat))
  			{
  				return (true);
  			}
  			// This is the contrast which appears on point 5 of the document.
  			return (Math.abs(contraststat) > limit);

  			// This is the contrast which appears on point 6 of the document.
  			// We need a function to convert the probability into the value of d*
  			// double val = - (this.nobservations / 2) * Math.Log (1 - Math.Pow (coefcorr, 2)) + 0.5;
  		}

  		/** 
  		 Returns if Values1 is correlated with Values2, being MinProb the % of 
  		 probability to achieve.
  		 
  		 @param Values1 Values1c= First Array of values with which calculate correlation.
  		 @param Values2 Values2 = Second Array of values with which calculate correlation.
  		 @param MinProb Percentage of probability with which the regressor must be selected.
  		 @return  if series are correlated 
  		*/
  		public static boolean Test_Correlation(double[] Values1, double[] Values2, double MinProb)
  		{
  			double coefcorr, contraststat, limit;
  			if (Values1.length != Values2.length)
  			{
  				String msg = String.format("Strings.Both_Vectors_must_have_the_same_length, Values1.length, Values2.length");
  				throw new IllegalArgumentException(msg);
  			}

  			limit = Stat.TStudent_quantil(MinProb / 100.0, Values1.length - 2, true);
  			coefcorr = CoefCorrelation(Values1, Values2);
  			contraststat = coefcorr / Math.sqrt((1.0 - Math.pow(coefcorr, 2.0)) / (Values1.length - 2.0));
  			if (Double.isNaN(contraststat))
  			{
  				return (true);
  			}
  			// This is the contrast which appears on point 5 of the document.
  			return (Math.abs(contraststat) > limit);
  		}

  		/** 
  		 Returns the value of the correlation contrast of Values1 being with Values2, being MinProb the % of 
  		 probability to achieve.
  		 
  		 @param Values1 Values1c= First Array of values with which calculate correlation.
  		 @param Values2 Values2 = Second Array of values with which calculate correlation.
  		 @param MinProb Percentage of probability with which the regressor must be selected.
  		 @return  result of contrast 
  		*/
  		public static double CorrelationContrast(double[] Values1, double[] Values2, double MinProb)
  		{
  			double coefcorr, contraststat, limit;
  			if (Values1.length != Values2.length)
  			{
  				String msg = String.format("Strings.Both_Vectors_must_have_the_same_length" + Values1.length + Values2.length);
  				throw new IllegalArgumentException(msg);
  			}

  			limit = Stat.TStudent_quantil(MinProb / 100.0, Values1.length - 2, true);
  			coefcorr = CoefCorrelation(Values1, Values2);
  			contraststat = coefcorr / Math.sqrt((1.0 - Math.pow(coefcorr, 2.0)) / (Values1.length - 2.0));
  			if (Double.isNaN(contraststat))
  			{
  				return (100.0);
  			}
  			// This is the contrast which appears on point 5 of the document.
  			return (Math.abs(contraststat) / limit);
  		}

  		/** Rounds trigonometric function to get the closest number to c.
  		 @param c c: Is the value to be rounded
  		 @return  result value 
  		*/
  		public static double RoundTrigFunc(double c)
  		{
  			double Val = c;

  			Val = (c < 0.00000001 && c > 0) ? (double) 0 : c;
  			Val = (c > 0.99999999 && c < 1) ? (double) 1 : c;
  			Val = (c > 0.49999999 && c < 0.5) ? (double) 0.5 : c;
  			Val = (c > 0.5 && c < 0.50000001) ? (double) 0.5 : c;

  			Val = (c > -0.00000001 && c < 0) ? (double) 0 : c;
  			Val = (c < -0.99999999 && c > -1) ? (double) - 1 : c;
  			Val = (c < -0.49999999 && c > -0.5) ? (double) - 0.5 : c;
  			Val = (c < -0.5 && c > -0.50000001) ? (double) - 0.5 : c;

  			return (Val);
  		}

  		/** 
  		 Returns the Correlation COeficient of array "Values1" and "Values2".
  		 
  		 @param Values1 Values1 = First Array of values from which to calculate the correlation coeficient.
  		 @param Values2 Values2 = Second Array of values from which to calculate the correlation coeficient.
  		 @return  calculated standard error 
  		*/
  		public static double StdError(double[] Values1, double[] Values2)
  		{
  			if (Values1.length != Values2.length)
  			{
  				throw new IllegalArgumentException("Strings.Both_vectors_must_have_the_same_number_of_values");
  			}
  			double Err = 0.0;

  			for (int i = 0;i < Values1.length;i++)
  			{
  				if (Values1[i] != 0.0)
  				{
  					Err += (Math.abs(Values1[i] - Values2[i])) / Values1[i];
  				}
  				else
  				{
  					Err += Math.abs(Values2[i]);
  				}
  			}
  			return Err / Values1.length;
  		}

  		//endregion Descriptive statistical function

  		//region Standardize Methods

  		/**  Bias standardizing of one value 
  		 @param X list of data for standardize 
  		 @param res bias 
  		 @param index index of value 
  		 @return  standardize value 
  		*/
  		public static double StandardizedRes(double[] X, double[] res, int index)
  		{
  			int n = X.length;
  			double Xmean = Stat.Mean(X, n, 0);
  			double hi = 1.0 / n + Math.pow((X[index] - Xmean), 2) / Stat.VarMuestral(X, n, 0);
  			double Sxy = 0.0;
  			for (int i = 0;i < n;i++)
  			{
  				Sxy += Math.pow(res[i], 2);
  			}
  			Sxy = Math.sqrt(Sxy / n - 2);

  			double SRi = res[index] / Sxy * Math.sqrt(1 - hi);
  			return SRi;
  		}

  		/**  Values standardizing 
  		 @param data list of data to standardize
  		 @param InflationFactor inflation factor 
  		 @return  standardized values 
  		*/
  		public static ArrayList<Double> StandarizeValues(ArrayList<Double> data, double InflationFactor) {
  			double mean = Stat.Mean(DotNetHelpers.DoubleLists.toArray(data), data.size(), 0);
  			double stddev = Math.sqrt(Stat.VarMuestral(DotNetHelpers.DoubleLists.toArray(data), data.size(), 0));
  			ArrayList<Double> v = new ArrayList<Double>();
  			for (int i = 0;i < data.size();i++) {
  				v.add((data.get(i) - mean / stddev) * InflationFactor);
  			}
  			return v;
  		}

  		//endregion

  		//region Robust statistics

  		//Obsolete methods. Better use Robust statistics region form Statistics class

  		/**  Trimmed Mean 
  		 @param values list of values 
  		 @param rec cut off in each boundaries 
  		 @return  trimmed mean 
  		*/
  		public static double TrimMean(double[] values, double rec)
  		{
  			int mitad = (int)((double)values.length / 2);
  			int nRec = (int)java.lang.Math.round((double)mitad * rec);
  			int nValues = values.length - nRec * 2;
  			if (nValues == 0)
  			{
  				return 0.0;
  			}
  			double sum = 0.0;
  			for (int i = nRec;i < nValues - nRec;i++)
  			{
  				sum += values[i];
  			}
  			return sum / nValues;
  		}

  		/**  Winsored Mean 
  		 @param values list of values 
  		 @param rec cut off in each boundaries 
  		 @return  winsored mean 
  		*/
  		public static double WinMean(double[] values, double rec)
  		{
  			int mitad = (int)((double)values.length / 2);
  			int nRec = (int)java.lang.Math.round((double)mitad * (1 - rec));
  			int nValues = values.length;
  			if (nValues == 0)
  			{
  				return 0.0;
  			}
  			double sum = 0.0;
  			for (int i = 0;i < nRec;i++)
  			{
  				sum += values[nRec];
  			}
  			for (int i = nRec;i < nValues - nRec;i++)
  			{
  				sum += values[i];
  			}
  			for (int i = nValues - nRec;i < nValues;i++)
  			{
  				sum += values[nValues - nRec - 1];
  			}
  			return sum / nValues;
  		}

  		//endregion

  }
  

    