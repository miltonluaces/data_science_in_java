package ClassicalStat;

import java.util.*;

//region Imports


//endregion


//DotNetHelpers_doc_comment_start  Hypothesis tests (parametric and non-parametric) 
//DotNetHelpers_doc_comment_body Results in p-values (probability of obtaining a result as extreme as the current, assuming null hypothesis is true, by chance). 
//DotNetHelpers_doc_comment_body Conventional interpretation: 
//DotNetHelpers_doc_comment_body p greater than 0.1      : Data consistent with null hypothesis
//DotNetHelpers_doc_comment_body p between 0.05 and 0.1  : Possibly null hypothesis does not hold. Mored data needed. Not conclusive.
//DotNetHelpers_doc_comment_body p between 0.01 and 0.05 : Some evidence against null hypothesis
//DotNetHelpers_doc_comment_body p less than 0.01        : Strong evidence against null hypothesis
//DotNetHelpers_doc_comment_body p less than 0.001       : Very strong evidence agains null hypothesis
//DotNetHelpers_doc_comment_end

public class HypoTests {

    //region Fields

    private Statistics stat;
    private NormalDist norm;
    private double significance;
    private int minCount;
    private HashMap<String, Double> AcCritValues;
    private double cutOff;


    //region Constants for normality test

    // Coefficients for P close to 0.5
    private double A0_p = 3.3871327179E+00, A1_p = 5.0434271938E+01, A2_p = 1.5929113202E+02, A3_p = 5.9109374720E+01, B1_p = 1.7895169469E+01, B2_p = 7.8757757664E+01, B3_p = 6.7187563600E+01;
    // Coefficients for P not close to 0, 0.5 or 1 (names changed to avoid conflict with Calculate)
    private double C0_p = 1.4234372777E+00, C1_p = 2.7568153900E+00, C2_p = 1.3067284816E+00, C3_p = 1.7023821103E-01, D1_p = 7.3700164250E-01, D2_p = 1.2021132975E-01;
    // Coefficients for P near 0 or 1.
    private double E0_p = 6.6579051150E+00, E1_p = 3.0812263860E+00, E2_p = 4.2868294337E-01, E3_p = 1.7337203997E-02, F1_p = 2.4197894225E-01, F2_p = 1.2258202635E-02;
    private double SPLIT1 = 0.425, SPLIT2 = 5.0, CONST1 = 0.180625, CONST2 = 1.6;

    // Constants & polynomial coefficients for alnorm(), slightly renamed to avoid conflicts.
    private double CON_a = 1.28, LTONE_a = 7.0, UTZERO_a = 18.66;
    private double P_a = 0.398942280444, Q_a = 0.39990348504, R_a = 0.398942280385, A1_a = 5.75885480458, A2_a = 2.62433121679, A3_a = 5.92885724438, B1_a = -29.8213557807, B2_a = 48.6959930692, C1_a = -3.8052E-8, C2_a = 3.98064794E-4, C3_a = -0.151679116635, C4_a = 4.8385912808, C5_a = 0.742380924027, C6_a = 3.99019417011, D1_a = 1.00000615302, D2_a = 1.98615381364, D3_a = 5.29330324926, D4_a = -15.1508972451, D5_a = 30.789933034;
    private double epsilon = 0.0001;
    //endregion


    //endregion

    //region Constructor

    //DotNetHelpers_doc_comment_start  Constructor 
    //DotNetHelpers_doc_comment_end
    public HypoTests(double significance, int minCount)
    {
        this.significance = significance;
        this.minCount = minCount;
        this.cutOff = 0;
        stat = new Statistics();
        norm = new NormalDist();
    }

    //endregion

    //region Properties

    //DotNetHelpers_doc_comment_start  test significance 
    //DotNetHelpers_doc_comment_end
    public final double getSignificance()
    {
        return significance;
    }

    //DotNetHelpers_doc_comment_start  minimum count for testing 
    //DotNetHelpers_doc_comment_end
    public final double getMinCount()
    {
        return minCount;
    }

    //DotNetHelpers_doc_comment_start  cut-off for trimming 
    //DotNetHelpers_doc_comment_end
    public final double getCutOff()
    {
        return cutOff;
    }
    public final void setCutOff(double value)
    {
        cutOff = value;
    }

    //endregion

    //region Parametric Tests

    //region Homogeneity of two sample distributions

    //DotNetHelpers_doc_comment_start  Parametric test for homogeneity of two samples 
    //DotNetHelpers_doc_comment_body @param S1 sample 1 
    //DotNetHelpers_doc_comment_body @param S2 sample 2 
    //DotNetHelpers_doc_comment_body @param pMean output mean comparison p-value 
    //DotNetHelpers_doc_comment_body @param pVar output variance comparison p-value 
    //DotNetHelpers_doc_comment_end
    public final void ParametricHomogeneity(ArrayList<Double> S1, List<Double> S2, DotNetHelpers.RefObject<Double> pMean, DotNetHelpers.RefObject<Double> pVar)
    {

        if (cutOff > 0)
        {
            S1 = stat.Trim((ArrayList<Double>)S1, 0.05, true);
            S2 = stat.Trim((ArrayList<Double>)S1, 0.05, true);
        }

        pVar.argValue = EqualVar(S1, S2);
        if(pVar.argValue < significance)
        {
            pMean.argValue = 0.0;
            return;
        }
        pMean.argValue = EqualMean(S1, S2);
    }

    //DotNetHelpers_doc_comment_start  Equal variance parametric test 
    //DotNetHelpers_doc_comment_body @param S1 sample 1 
    //DotNetHelpers_doc_comment_body @param S2 sample 2 
    //DotNetHelpers_doc_comment_body @return  p-value 
    //DotNetHelpers_doc_comment_end
    public final double EqualVar(List<Double> S1, List<Double> S2)
    {
        double v1 = stat.Variance(S1);
        double v2 = stat.Variance(S2);
        if(S1.size() < minCount || S2.size() < minCount || Math.abs(v1-v2) < epsilon)
        {
            return 1.0;
        }
        int df1 = S1.size()-1;
        int df2 = S2.size()-1;
        double f;
        double p = -1;
        if(v1 > v2)
        {
            f = v1/v2;
            if(Math.abs(f) > 10)
            {
                return 0;
            }
            p = Stat.FDistribution(f, df1, df2);
        }
        else
        {
            f = v2/v1;
            if(Math.abs(f) > 10)
            {
                return 0;
            }
            p = Stat.FDistribution(f, df2, df1);
        }
        p = 2 * p;
        if(p > 1)
        {
            return 1;
        }
        return p;
    }

   //DotNetHelpers_doc_comment_start  Equal mean parametric test 
    //DotNetHelpers_doc_comment_body @param S1 sample 1 
    //DotNetHelpers_doc_comment_body @param S2 sample 2 
    //DotNetHelpers_doc_comment_body @return  p-value 
   //DotNetHelpers_doc_comment_end
    public final double EqualMean(List<Double> S1, List<Double> S2)
    {
        double m1 = stat.Mean(S1);
        double m2 = stat.Mean(S2);
        double n1 = S1.size();
        double n2 = S2.size();

        if(n1 < minCount || n2 < minCount)
        {
            double pn = -1;
            if(m1 > m2)
            {
                pn = 1- norm.pNorm(1 - m2/m1, 0, m2/m1);
            }
            else
            {
                pn = 1- norm.pNorm(1 - m1/m2, 0, m1/m2);
            }
            if(pn > 1)
            {
                return 1;
            }
            return pn;
        }

        ArrayList<Double> S = new ArrayList<Double>(S1);
        S.addAll(S2);
        double var = stat.Variance(S);

        int df = (int)(n1+n2-2);
        double t = (m1 - m2) / Math.sqrt(var * (1/n1 + 1/n2));
        if(t == 0)
        {
            return 1;
        }
        //double p = 2 * Stat.TStudent_acum (t, 0, 0, df, false);
        //double p = Stat.TStudent_acum(t, 0, 0, df, false);
        double p = stat.pt(t, df);

        if(p > 1)
        {
            return 1;
        }
        return p;
    }

    //endregion

    //region Homoscedasticity

    //region Bartlett

    public final double BartlettTest(List<List<Double>> groups)
    {
        double x2 = GetX2(groups);
        int df = groups.size() - 1;
        return 1-Stat.chi_square_p(x2, df);
    }

    private double GetX2(List<List<Double>> groups)
    {
        int kInt = groups.size();
        double k = (double)kInt;
        double[] n = new double[kInt];
        double N = 0;
        double[] v = new double[kInt];
        for (int i = 0; i < groups.size(); i++)
        {
            n[i] = groups.get(i).size();
            N += n[i];
            v[i] = stat.Variance(groups.get(i));
        }

        //Sp = 1/(N-k) Sum (ni-1)Vi
        double Sp = 0;
        for (int i = 0; i < k; i++)
        {
            Sp += (n[i] - 1) * v[i];
        }
        Sp *= 1 / (N - k);

        //Sl = Sum (ni-1)ln(Vi)
        double Sl = 0;
        for (int i = 0; i < k; i++)
        {
            Sl += (n[i] - 1) * Math.log(v[i]);
        }

        //Sq = Sum 1/(ni-1)
        double Sq = 0;
        for (int i = 0; i < k; i++)
        {
            Sq += (1 / (n[i] - 1)) - 1/(N-k);
        }

        double num = (N - k) * Math.log(Sp) - Sl;
        double den = 1 + 1 / (3 * (k - 1)) * Sq;

        return num / den;
    }

    //endregion

    //region Kruskal Wallis

    public final double KruskalWallis(List<List<Double>> groups, boolean tiedCorrection)
    {
        double KW = GetKW(groups, tiedCorrection);
        int df = groups.size() - 1;
        return 1-Stat.chi_square_p(KW, df);
    }

    private double GetKW(List<List<Double>> groups, boolean tiedCorrection)
    {
        // data groups, number of elements ni
        double g = groups.size();
        ArrayList<Double> n = new ArrayList<Double>();
        ArrayList<ArrayList<Datum>> dGroups = new ArrayList<ArrayList<Datum>>();
        ArrayList<Datum> dGroup;
        Datum datum;
        for (int i = 0; i < groups.size(); i++)
        {
            dGroup = new ArrayList<Datum>();
            n.add((double) groups.get(i).size());
            for (int j = 0; j < groups.get(i).size(); j++)
            {
                datum = new Datum(groups.get(i).get(j), i);
                dGroup.add(datum);
            }
            dGroups.add(dGroup);
        }

        //general group data : ranking
        ArrayList<Datum> data = new ArrayList<Datum>();
        for (ArrayList<Datum> dg : dGroups)
        {
            data.addAll(dg);
        }
        Collections.sort(data);

        double N = data.size();
        double rank = 1;
        data.get((int)N-1).rank = 1;
        for (int k = (int)N-1; k > 0; k--)
        {
            if (data.get(k).value > data.get(k - 1).value)
            {
                rank++;
            }
            data.get(k).rank = rank;
        }

        //mean ranking foreach group
        ArrayList<Double> r = new ArrayList<Double>();
        double rSum = 0;
        for (ArrayList<Datum> dg : dGroups)
        {
            for (Datum d : dg)
            {
                rSum += d.rank;
            }
            r.add(rSum / (double)dg.size());
        }

        //mean ranking of all:
        double R = (N + 1) / 2;

        //KW statistic
        double num = 0;
        for (int i = 0; i < g; i++)
        {
            num += n.get(i) * Math.pow(r.get(i) - R, 2);
        }

        double den = 0;
        for (int i = 0; i < g; i++)
        {
            for (int j = 0; j < n.get(i); j++)
            {
                den += n.get(i) * Math.pow(dGroups.get(i).get(j).rank - R, 2);
            }
        }

        double KW = (N - 1) * num / den;
        if (tiedCorrection)
        {
            KW = KWTiedCorrection(dGroups, KW, N);
        }
        return KW;
    }

    private double KWTiedCorrection(ArrayList<ArrayList<Datum>> dGroups, double kw, double N)
    {
        double[] t = new double[dGroups.size()];
        for (int i = 0; i < dGroups.size(); i++)
        {
            t[i] = NTied(dGroups.get(i));
        }

        double C = 0;
        for (int i = 0; i < dGroups.size(); i++)
        {
            C += Math.pow(t[i], 3) - t[i];
        }
        C /= (Math.pow(N, 3) - N);
        return kw / (1 - C);
    }

    private double NTied(ArrayList<Datum> dGroup)
    {
        double t = 0;
        Collections.sort(dGroup);
        for (int j = 1; j < dGroup.size();j++)
        {
            if (dGroup.get(j).rank == dGroup.get(j - 1).rank)
            {
                t++;
            }
        }
        return t;
    }

    //endregion

    //endregion

    //endregion

    //region Non-Parametric Tests

    //region Homogeneity of two sample distributions

    //region Mann-Withney ranks test (means)

    //DotNetHelpers_doc_comment_start  Non-parametric ranges test for homogeneity of samples 
    //DotNetHelpers_doc_comment_body @param S1 sample 1 
    //DotNetHelpers_doc_comment_body @param S2 sample 2 
    //DotNetHelpers_doc_comment_body @param twoSided true if it is a two-sided test, false if it is one-sided 
    //DotNetHelpers_doc_comment_body @return  p-value 
    //DotNetHelpers_doc_comment_end
    public final double MannWhitney(List<Double> S1, List<Double> S2, boolean twoSided)
    {
        ArrayList<Datum> data = new ArrayList<Datum>();
        for(double val : S1)
        {
            data.add(new Datum(val, 1));
        }
        for(double val : S2)
        {
            data.add(new Datum(val, 2));
        }
        double n1 = (double)S1.size();
        double n2 = (double)S2.size();
        return MannWhitney(data, n1, n2, twoSided);
    }

    //DotNetHelpers_doc_comment_start  Non-parametric ranges test for homogeneity of a partition of sample  
    //DotNetHelpers_doc_comment_body @param sample
    //DotNetHelpers_doc_comment_body @param index first index of second part of the sample 
    //DotNetHelpers_doc_comment_body @param twoSided if the hypothesis test is two-sided or one-sided 
    //DotNetHelpers_doc_comment_body @return  p-value 
    //DotNetHelpers_doc_comment_end
    public final double MannWhitney(ArrayList<Double> sample, int index, boolean twoSided)
    {
        ArrayList<Datum> data = new ArrayList<Datum>();
        for(int i=0;i<index;i++)
        {
            data.add(new Datum(sample.get(i), 1));
        }
        for(int i=index;i<sample.size();i++)
        {
            data.add(new Datum(sample.get(i), 2));
        }
        double n1 = (double)index;
        double n2 = (double)sample.size() - index;
        return MannWhitney(data, n1, n2, twoSided);
    }

    private double MannWhitney(ArrayList<Datum> data, double n1, double n2, boolean twoSided)
    {
        if(n1 == 0 || n2 == 0)
        {
            return 1;
        }
        //calculate ranges and R1, with tie-breaks
        Collections.sort(data);
        double rank = 1;
        double sumTie;
        double promTie;
        double R1 = 0;
        int tie = -1;
        for(int i=0;i<data.size();i++)
        {
            data.get(i).rank = rank++;
            if(i > 0 && tie == -1 && data.get(i).value == data.get(i-1).value)
            {
                tie = i-1;
            }
            if(i == data.size()-1 && tie != -1)
            {
                sumTie = 0;
                for(int j=tie;j<=i;j++)
                {
                    sumTie += data.get(j).rank;
                }
                promTie = sumTie / (i-tie+1);
                for(int j=tie;j<=i;j++)
                {
                    data.get(j).rank = promTie;
                }
                tie = -1;
            }
            if(tie != -1 && data.get(i).value != data.get(i-1).value)
            {
                sumTie = 0;
                for(int j=tie;j<i;j++)
                {
                    sumTie += data.get(j).rank;
                }
                promTie = sumTie / (i-tie);
                for(int j=tie;j<i;j++)
                {
                    data.get(j).rank = promTie;
                }
                tie = -1;
            }
        }

        //calculate R1, U statistic, Expected value and variance
        for(Datum d : data)
        {
            if(d.sample == 1)
            {
                R1 += d.rank;
            }
        }
        double U = n1 * n2 + (n1 * (n1 + 1))/ 2 - R1;
        double EU = (n1 * n2) / 2.0;
        double VarU = (n1 * n2 * (n1 + n2 + 1))/12.0;

        //calculate z statistic and p-value
        double p = 1- norm.pNorm(U, EU, Math.sqrt(VarU));
        if(twoSided)
        {
            p = 2 * p;
        }
        if(p > 1)
        {
            p = 1;
        }
        return p;
    }

    //endregion

    //region Wald-Wolfowitz runs test (means, variances, simmetries)

    //DotNetHelpers_doc_comment_start  Non-parametric runs test for homogeneity of samples 
    //DotNetHelpers_doc_comment_body @param S1 sample 1 
    //DotNetHelpers_doc_comment_body @param S2 sample 2 
    //DotNetHelpers_doc_comment_body @param twoSided true if it is a two-sided test, false if it is one-sided 
    //DotNetHelpers_doc_comment_body @return  p-value 
    //DotNetHelpers_doc_comment_end
    public final double WaldWolfowitz(List<Double> S1, List<Double> S2, boolean twoSided)
    {

        //load data
        ArrayList<Datum> data = new ArrayList<Datum>();
        for(double val : S1)
        {
            data.add(new Datum(val, 1));
        }
        for(double val : S2)
        {
            data.add(new Datum(val, 2));
        }
        double n1 = (double)S1.size();
        double n2 = (double)S2.size();

        //calculate runs
        Collections.sort(data);
        int R = 0;
        int sAnt = -1;
        for(Datum d : data)
        {
            if(d.sample != sAnt)
            {
                R++;
            }
            sAnt = d.sample;
        }

        //calculate mean and variance
        double ER = ((2 * n1 * n2) / (n1 + n2)) + 1;
        double VarR = (2 * n1 * n2 * (2 * n1 *n2 - n1 - n2)) / (Math.pow((n1 + n2), 2) * (n1 + n2 - 1));
        double p = norm.pNorm(R, ER, VarR);
        if(twoSided)
        {
            p = 2 * p;
        }
        if(p > 1)
        {
            p = 1;
        }
        return p;
    }

    //endregion

    //endregion

    //region Normality of a sample

    //region Shapiro-Wilk Test

    //region Public Method

    //DotNetHelpers_doc_comment_start  Shapiro-Wilk. Returns an array of size two if everything goes well, or null if there's an error. 
    //DotNetHelpers_doc_comment_body @param sample sample to test normality 
    //DotNetHelpers_doc_comment_body @return  p value 
    //DotNetHelpers_doc_comment_end
    public final double ShapiroWilk(List<Double> sampleIList)
    {
        if (sampleIList.size() < 3)
        {
            return -1;
        }
        ArrayList<Double> sampleList = new ArrayList<Double>(sampleIList);
        Collections.sort(sampleList);
        double[] sample = new double[sampleList.size()+1];
        for(int i=1;i<sample.length;i++)
        {
            sample[i] = sampleList.get(i-1);
        }

        //constants
        double[] C1 = { Double.NaN, 0.0E0, 0.221157E0, -0.147981E0, -0.207119E1, 0.4434685E1, -0.2706056E1 };
        double[] C2 = { Double.NaN, 0.0E0, 0.42981E-1, -0.293762E0, -0.1752461E1, 0.5682633E1, -0.3582633E1 };
        double[] C3 = { Double.NaN, 0.5440E0, -0.39978E0, 0.25054E-1, -0.6714E-3 };
        double[] C4 = { Double.NaN, 0.13822E1, -0.77857E0, 0.62767E-1, -0.20322E-2 };
        double[] C5 = { Double.NaN, -0.15861E1, -0.31082E0, -0.83751E-1, 0.38915E-2 };
        double[] C6 = { Double.NaN, -0.4803E0, -0.82676E-1, 0.30302E-2 };
        double[] C7 = { Double.NaN, 0.164E0, 0.533E0 };
        double[] C8 = { Double.NaN, 0.1736E0, 0.315E0 };
        double[] C9 = { Double.NaN, 0.256E0, -0.635E-2 };
        double[] G = { Double.NaN, -0.2273E1, 0.459E0 };
        double z90 = 0.12816E1;
        double z95 = 0.16449E1;
        double z99 = 0.23263E1;
        double zm = 0.17509E1;
        double zss = 0.56268E0;
        double bf1 = 0.8378E0;
        double xx90 = 0.556E0;
        double xx95 = 0.622E0;
        double sqrth = 0.70711E0;
        double th = 0.375E0;
        double epsilon = 1E-19;
        double pi6 = 0.1909859E1;
        double stqr = 0.1047198E1;
        boolean upper = true;

        int n1 = sample.length-1;
        int n2 = (sample.length-1)/2;
        double[] a = new double[sample.length];
        double[] w = new double[1];
        double[] pw = new double[1];
        int result = -1;
        boolean init = false;
        int n = sample.length-1;
        pw[0] = 1.0;
        if(w[0] >= 0.0)
        {
            w[0] = 1.0;
        }
        double an = n;
        result = 3;
        int nn2 = n/2;
        if(n2 < nn2)
        {
            return pw[0];
        }
        result = 1;
        if(n < 3)
        {
            return pw[0];
        }

        // If init is false, calculates coefficients for the test
        if(!init)
        {
            if(n == 3)
            {
                a[1] = sqrth;
            }
            else
            {
                double an25 = an + 0.25;
                double summ2 = 0.0;
                for(int i = 1;i <= n2;++i)
                {
                    a[i] = ppnd((i - th) / an25);
                    summ2 += a[i] * a[i];
                }
                summ2 *= 2.0;
                double ssumm2 = Math.sqrt(summ2);
                double rsn = 1.0 / Math.sqrt(an);
                double a1 = poly(C1, 6, rsn) - a[1] / ssumm2;

                // Normalize coefficients

                int i1;
                double fac;
                if(n > 5)
                {
                    i1 = 3;
                    double a2 = -a[2] / ssumm2 + poly(C2, 6, rsn);
                    fac = Math.sqrt((summ2 - 2.0 * a[1] * a[1] - 2.0 * a[2] * a[2]) / (1.0 - 2.0 * a1 * a1 - 2.0 * a2 * a2));
                    a[1] = a1;
                    a[2] = a2;
                }
                else
                {
                    i1 = 2;
                    fac = Math.sqrt((summ2 - 2.0 * a[1] * a[1]) / (1.0 - 2.0 * a1 * a1));
                    a[1] = a1;
                }
                for(int i = i1;i <= nn2;++i)
                {
                    a[i] = -a[i] / fac;
                }
            }
            init = true;
        }
        if(n1 < 3)
        {
            return pw[0];
        }
        int ncens = n - n1;
        result = 4;
        if(ncens < 0 || (ncens > 0 && n < 20))
        {
            return pw[0];
        }
        result = 5;
        double delta = ncens / an;
        if(delta > 0.8)
        {
            return pw[0];
        }

        // If W input as negative, calculate significance level of -W

        double w1, xx;
        if(w[0] < 0.0)
        {
            w1 = 1.0 + w[0];
            result = 0;
        }
        else
        {

            // Check for zero range
            result = 6;
            double range = sample[n1] - sample[1];
            if(range < epsilon)
            {
                return pw[0];
            }

            // Check for correct sort order on range - scaled X
            result = 7;
            xx = sample[1] / range;
            double sx = xx;
            double sa = -a[1];
            int j = n - 1;
            for(int i=2;i<=n1;++i)
            {
                double xi = sample[i] / range;
                // IF (XX-XI .GT. epsilon) PRINT *,' ANYTHING'
                sx += xi;
                if(i != j)
                {
                    sa += sign(1, i - j) * a[Math.min(i, j)];
                }
                xx = xi;
                --j;
            }
            result = 0;
            if(n > 5000)
            {
                result = 2;
            }

            // Calculate W statistic as squared correlation between data and coefficients
            sa /= n1;
            sx /= n1;
            double ssa = 0.0;
            double ssx = 0.0;
            double sax = 0.0;
            j = n;
            double asa;
            for(int i = 1;i <= n1;++i)
            {
                if(i != j)
                {
                    asa = sign(1, i - j) * a[Math.min(i, j)] - sa;
                }
                else
                {
                    asa = -sa;
                }
                double xsx = sample[i] / range - sx;
                ssa += asa * asa;
                ssx += xsx * xsx;
                sax += asa * xsx;
                --j;
            }

            // W1 equals (1-W) calculated to avoid excessive rounding error for W very near 1 (a potential problem in very large samples)
            double ssassx = Math.sqrt(ssa * ssx);
            w1 = (ssassx - sax) * (ssassx + sax) / (ssa * ssx);
        }
        w[0] = 1.0 - w1;

        // Calculate significance level for W (exact for N=3)
        if(n == 3)
        {
            pw[0] = pi6 * (Math.asin(Math.sqrt(w[0])) - stqr);
            return pw[0];
        }
        double y = Math.log(w1);
        xx = Math.log(an);
        double m = 0.0;
        double s = 1.0;
        if(n <= 11)
        {
            double gamma = poly(G, 2, an);
            if(y >= gamma)
            {
                pw[0] = epsilon;
                return pw[0];
            }
            y = -Math.log(gamma - y);
            m = poly(C3, 4, an);
            s = Math.exp(poly(C4, 4, an));
        }
        else
        {
            m = poly(C5, 4, xx);
            s = Math.exp(poly(C6, 3, xx));
        }
        if(ncens > 0)
        {

            // Censoring by proportion NCENS/N. Calculate mean and sd of normal equivalent deviate of W.
            double ld = -Math.log(delta);
            double bf = 1.0 + xx * bf1;
            double z90f = z90 + bf * Math.pow(poly(C7, 2, Math.pow(xx90, xx)), ld);
            double z95f = z95 + bf * Math.pow(poly(C8, 2, Math.pow(xx95, xx)), ld);
            double z99f = z99 + bf * Math.pow(poly(C9, 2, xx), ld);

            // Regress Z90F,...,Z99F on normal deviates z90,...,z99 to get pseudo-mean and pseudo-sd of z as the slope and intercept
            double zfm = (z90f + z95f + z99f) / 3.0;
            double zsd = (z90 * (z90f - zfm) + z95 * (z95f - zfm) + z99 * (z99f - zfm)) / zss;
            double zbar = zfm - zsd * zm;
            m += zbar * s;
            s *= zsd;
        }
        pw[0] = alnorm((y - m) / s, upper);
        if(result != 0 && result != 2)
        {
            return (double)result;
        }
        return pw[0];
    }

    //endregion

    //region Private Methods

    //Returns an int with abs value of x and sign of y 
    private int sign(int value, int sign)
    {
        int result = Math.abs(value);
        if(sign < 0)
        {
            return -result;
        }
        return result;
    }

    //Produces the normal deviate Z corresponding to a given lower tail area of P
    private double ppnd(double p)
    {
        double q = p - 0.5;
        double r;
        if(Math.abs(q) <= SPLIT1)
        {
            r = CONST1 - q * q;
            return q * (((A3_p * r + A2_p) * r + A1_p) * r + A0_p) / (((B3_p * r + B2_p) * r + B1_p) * r + 1.0);
        }
        else
        {
            if(q < 0.0)
            {
                r = p;
            }
            else
            {
                r = 1.0 - p;
            }
            if(r <= 0.0)
            {
                return 0.0;
            }

            r = Math.sqrt(-Math.log(r));
            double normal_dev;
            if(r <= SPLIT2)
            {
                r -= CONST2;
                normal_dev = (((C3_p * r + C2_p) * r + C1_p) * r + C0_p) / ((D2_p * r + D1_p) * r + 1.0);
            }
            else
            {
                r -= SPLIT2;
                normal_dev = (((E3_p * r + E2_p) * r + E1_p) * r + E0_p) / ((F2_p * r + F1_p) * r + 1.0);
            }
            if(q < 0.0)
            {
                normal_dev = -normal_dev;
            }
            return normal_dev;
        }
    }

    //Evaluates polynom of order nord-1 with array of coefficients c. Zero order coefficient is c[1]
    private double poly(double[] c, int nord, double x)
    {
        double poly = c[1];
        if(nord == 1)
        {
            return poly;
        }
        double p = x * c[nord];
        if(nord != 2)
        {
            int n2 = nord - 2;
            int j = n2 + 1;
            for(int i = 1;i <= n2;++i)
            {
                p = (p + c[j]) * x;
                --j;
            }
        }
        poly += p;
        return poly;
    }

    //Evaluates the tail area of the standardised normal curve from x to infinity if upper is true or from minus infinity to sample if upper is false.
    private double alnorm(double x, boolean upper)
    {
        boolean up = upper;
        double z = x;
        if(z < 0.0)
        {
            up = !up;
            z = -z;
        }
        double fn_val;
        if(z > LTONE_a && (!up || z > UTZERO_a))
        {
            fn_val = 0.0;
        }
        else
        {
            double y = 0.5 * z * z;
            if(z <= CON_a)
            {
                fn_val = 0.5 - z * (P_a - Q_a * y / (y + A1_a + B1_a / (y + A2_a + B2_a / (y + A3_a))));
            }
            else
            {
                fn_val = R_a * Math.exp(-y) / (z + C1_a + D1_a / (z + C2_a + D2_a / (z + C3_a + D3_a / (z + C4_a + D4_a / (z + C5_a + D5_a / (z + C6_a))))));
            }
        }
        if(!up)
        {
            fn_val = 1.0 - fn_val;
        }
        return fn_val;
    }

    //endregion

    //endregion

    //endregion

    //endregion

    //region ANOVA for Nested Models

    //DotNetHelpers_doc_comment_start  ANOVA for nested models test  
    //DotNetHelpers_doc_comment_body @param data real time series data 
    //DotNetHelpers_doc_comment_body @param redModel reduced model 
    //DotNetHelpers_doc_comment_body @param compModel complete model 
    //DotNetHelpers_doc_comment_body @param pr reduced model number of parameters 
    //DotNetHelpers_doc_comment_body @param pc complete model number of parameters
    //DotNetHelpers_doc_comment_body @param alpha significance level 
    //DotNetHelpers_doc_comment_body @return  if the new parameter(s) (is)are relevant 
    //DotNetHelpers_doc_comment_end
    public final boolean AnovaNestedModels(List<Double> data, List<Double> redModel, ArrayList<Double> compModel, int pr, int pc, double alpha)
    {
        if(data.size() != redModel.size() || data.size() != compModel.size())
        {
            throw new RuntimeException("Strings.Data_must_have_the_same_size_as_model");
        }
        double sser = SSE(data, redModel);
        double ssec = SSE(data, compModel);
        int n = data.size();
        if(sser == 0)
        {
            return false;
        }
        if(ssec == 0)
        {
            return true;
        }

        double F = ((sser - ssec) / (pc-pr)) / (ssec / (n - (pc+1)));
        if(F < 0)
        {
            return false;
        }
        if(ssec < sser && ssec == 0)
        {
            return true;
        }
        int v1 = pc-pr;
        int v2 = n - (pc-1);

        return Stat.Test_F(F, v1, v2) > 1-alpha;
    }

    private double SSE(List<Double> data, List<Double> model)
    {
        if(data.size() != model.size())
        {
            throw new RuntimeException("Strings.Data_must_have_the_same_size_as_model");
        }

        double sse = 0.0;
        for(int i=0;i<data.size();i++)
        {
            sse += Math.pow(data.get(i) - model.get(i), 2);
        }
        return sse;
    }

    //endregion


		//region ACF tests
		//region Durbin-Watson
		/**  Parametric autocorrelation test. Assumption: normality of residuals (not robust) 
		 d = 2: no correlation. d less than 2 : possitive autocorrelation. if d less than 1 great positive autocorrelation. idem for negative
		 @param serie time series 
		 @param model model for the time series (if null, test (not recommended if its not normal) on the time series 
		 @return  d value */

public double GetDWDValue(List<Double> serie, List<Double> model) {

	ArrayList<Double> res = new ArrayList<Double>(serie);
	if(model != null) {
		for(int i=1;i<res.size();i++) {
			//TODO: see what original code was here
		}
	}
	return -1;
}

public double GetDWCritValue(int n, int p, Boolean lower, Boolean signFive) {
if (AcCritValues == null) {	InitAutocorrelationTest();	}

String key = "";

//number of cases
if (n <= 15)
{
	key += "15-";
}
else if (n <= 20)
{
	key += "20-";
}
else if (n <= 25)
{
	key += "25-";
}
else if (n <= 30)
{
	key += "30-";
}
else if (n <= 40)
{
	key += "40-";
}
else if (n <= 50)
{
	key += "50-";
}
else if (n <= 60)
{
	key += "60-";
}
else if (n <= 80)
{
	key += "80-";
}
else
{
	key += "100-";
}

//number of parameters
if (p <= 2)
{
	key += "2-";
}
else if (p <= 3)
{
	key += "3-";
}
else if (p <= 4)
{
	key += "4-";
}
else if (p <= 5)
{
	key += "5-";
}
else
{
	key += "6-";
}

//significance 0.05 or 0.01
if (signFive)
{
	key += "0.05-";
}
else
{
	key += "0.01";
}

//lower or upper critical value
if (lower)
{
	key += "l";
}
else
{
	key += "u";
}

if (!AcCritValues.containsKey(key))
{
	throw new RuntimeException("Strings.Key_0_not_found, key");
}
return AcCritValues.get(key);
}

/**  Test Durbin-Watson of autocorrelation 
@param serie time series 
@param model model for the time series (if null, autocorrelation of time series, not recommended) 
@param lag lag for autocorrelation 
@param p number of model parameters 
@param signFive significance level = 0.05 (if false, 0.01) 
@return  result: 0: no correlation, -1: negative correlation. +1 positive correlation 
*/
//C# TO JAVA CONVERTER TODO TASK: The following line could not be converted:
public int IsAutocorrelated(List<Double> serie, List<Double> model, int p, Boolean signFive)
{
double d = GetDWDValue(serie, model);
int n = serie.size();
double cvl = GetDWCritValue(n, p, true, signFive);
double cvu = GetDWCritValue(n, p, false, signFive);

if (d < cvl)
{
	return +1;
}
else if (d > cvu)
{
	return -1;
}
else
{
	return 0;
}
}

/**  Initiallize autocorrelation test 
*/
//C# TO JAVA CONVERTER TODO TASK: The following line could not be converted:
public void InitAutocorrelationTest()
{
AcCritValues = new HashMap<String, Double>();
AcCritValues.put("15-2-0.05-l", 1.08);
AcCritValues.put("15-2-0.05-u", 1.36);
AcCritValues.put("15-3-0.05-l", 0.95);
AcCritValues.put("15-3-0.05-u", 1.54);
AcCritValues.put("15-4-0.05-l", 0.82);
AcCritValues.put("15-4-0.05-u", 1.75);
AcCritValues.put("15-5-0.05-l", 0.69);
AcCritValues.put("15-5-0.05-u", 1.97);
AcCritValues.put("15-6-0.05-l", 0.56);
AcCritValues.put("15-6-0.05-u", 2.21);
AcCritValues.put("15-2-0.01-l", 0.81);
AcCritValues.put("15-2-0.01-u", 1.07);
AcCritValues.put("15-3-0.01-l", 0.7);
AcCritValues.put("15-3-0.01-u", 1.25);
AcCritValues.put("15-4-0.01-l", 0.59);
AcCritValues.put("15-4-0.01-u", 1.46);
AcCritValues.put("15-5-0.01-l", 0.49);
AcCritValues.put("15-5-0.01-u", 1.7);
AcCritValues.put("15-6-0.01-l", 0.39);
AcCritValues.put("15-6-0.01-u", 1.96);

AcCritValues.put("20-2-0.05-l", 1.2);
AcCritValues.put("20-2-0.05-u", 1.41);
AcCritValues.put("20-3-0.05-l", 1.1);
AcCritValues.put("20-3-0.05-u", 1.54);
AcCritValues.put("20-4-0.05-l", 1.0);
AcCritValues.put("20-4-0.05-u", 1.68);
AcCritValues.put("20-5-0.05-l", 0.9);
AcCritValues.put("20-5-0.05-u", 1.83);
AcCritValues.put("20-6-0.05-l", 0.79);
AcCritValues.put("20-6-0.05-u", 1.99);
AcCritValues.put("20-2-0.01-l", 0.95);
AcCritValues.put("20-2-0.01-u", 1.15);
AcCritValues.put("20-3-0.01-l", 0.86);
AcCritValues.put("20-3-0.01-u", 1.27);
AcCritValues.put("20-4-0.01-l", 0.77);
AcCritValues.put("20-4-0.01-u", 1.41);
AcCritValues.put("20-5-0.01-l", 0.68);
AcCritValues.put("20-5-0.01-u", 1.57);
AcCritValues.put("20-6-0.01-l", 0.6);
AcCritValues.put("20-6-0.01-u", 1.74);

AcCritValues.put("25-2-0.05-l", 1.29);
AcCritValues.put("25-2-0.05-u", 1.45);
AcCritValues.put("25-3-0.05-l", 1.21);
AcCritValues.put("25-3-0.05-u", 1.55);
AcCritValues.put("25-4-0.05-l", 1.12);
AcCritValues.put("25-4-0.05-u", 1.66);
AcCritValues.put("25-5-0.05-l", 1.04);
AcCritValues.put("25-5-0.05-u", 1.77);
AcCritValues.put("25-6-0.05-l", 0.95);
AcCritValues.put("25-6-0.05-u", 1.89);
AcCritValues.put("25-2-0.01-l", 1.05);
AcCritValues.put("25-2-0.01-u", 1.21);
AcCritValues.put("25-3-0.01-l", 0.98);
AcCritValues.put("25-3-0.01-u", 1.3);
AcCritValues.put("25-4-0.01-l", 0.9);
AcCritValues.put("25-4-0.01-u", 1.41);
AcCritValues.put("25-5-0.01-l", 0.83);
AcCritValues.put("25-5-0.01-u", 1.52);
AcCritValues.put("25-6-0.01-l", 0.75);
AcCritValues.put("25-6-0.01-u", 1.65);

AcCritValues.put("30-2-0.05-l", 1.35);
AcCritValues.put("30-2-0.05-u", 1.49);
AcCritValues.put("30-3-0.05-l", 1.28);
AcCritValues.put("30-3-0.05-u", 1.57);
AcCritValues.put("30-4-0.05-l", 1.21);
AcCritValues.put("30-4-0.05-u", 1.65);
AcCritValues.put("30-5-0.05-l", 1.14);
AcCritValues.put("30-5-0.05-u", 1.74);
AcCritValues.put("30-6-0.05-l", 1.07);
AcCritValues.put("30-6-0.05-u", 1.83);
AcCritValues.put("30-2-0.01-l", 1.13);
AcCritValues.put("30-2-0.01-u", 1.26);
AcCritValues.put("30-3-0.01-l", 1.07);
AcCritValues.put("30-3-0.01-u", 1.34);
AcCritValues.put("30-4-0.01-l", 1.01);
AcCritValues.put("30-4-0.01-u", 1.42);
AcCritValues.put("30-5-0.01-l", 1.94);
AcCritValues.put("30-5-0.01-u", 1.51);
AcCritValues.put("30-6-0.01-l", 1.88);
AcCritValues.put("30-6-0.01-u", 1.61);

AcCritValues.put("40-2-0.05-l", 1.44);
AcCritValues.put("40-2-0.05-u", 1.54);
AcCritValues.put("40-3-0.05-l", 1.39);
AcCritValues.put("40-3-0.05-u", 1.60);
AcCritValues.put("40-4-0.05-l", 1.34);
AcCritValues.put("40-4-0.05-u", 1.66);
AcCritValues.put("40-5-0.05-l", 1.39);
AcCritValues.put("40-5-0.05-u", 1.72);
AcCritValues.put("40-6-0.05-l", 1.23);
AcCritValues.put("40-6-0.05-u", 1.79);
AcCritValues.put("40-2-0.01-l", 1.25);
AcCritValues.put("40-2-0.01-u", 1.34);
AcCritValues.put("40-3-0.01-l", 1.20);
AcCritValues.put("40-3-0.01-u", 1.40);
AcCritValues.put("40-4-0.01-l", 1.15);
AcCritValues.put("40-4-0.01-u", 1.46);
AcCritValues.put("40-5-0.01-l", 1.40);
AcCritValues.put("40-5-0.01-u", 1.52);
AcCritValues.put("40-6-0.01-l", 1.05);
AcCritValues.put("40-6-0.01-u", 1.58);

AcCritValues.put("50-2-0.05-l", 1.50);
AcCritValues.put("50-2-0.05-u", 1.59);
AcCritValues.put("50-3-0.05-l", 1.46);
AcCritValues.put("50-3-0.05-u", 1.63);
AcCritValues.put("50-4-0.05-l", 1.42);
AcCritValues.put("50-4-0.05-u", 1.67);
AcCritValues.put("50-5-0.05-l", 1.38);
AcCritValues.put("50-5-0.05-u", 1.72);
AcCritValues.put("50-6-0.05-l", 1.34);
AcCritValues.put("50-6-0.05-u", 1.77);
AcCritValues.put("50-2-0.01-l", 1.32);
AcCritValues.put("50-2-0.01-u", 1.40);
AcCritValues.put("50-3-0.01-l", 1.28);
AcCritValues.put("50-3-0.01-u", 1.45);
AcCritValues.put("50-4-0.01-l", 1.24);
AcCritValues.put("50-4-0.01-u", 1.49);
AcCritValues.put("50-5-0.01-l", 1.20);
AcCritValues.put("50-5-0.01-u", 1.54);
AcCritValues.put("50-6-0.01-l", 1.16);
AcCritValues.put("50-6-0.01-u", 1.59);

AcCritValues.put("60-2-0.05-l", 1.55);
AcCritValues.put("60-2-0.05-u", 1.62);
AcCritValues.put("60-3-0.05-l", 1.51);
AcCritValues.put("60-3-0.05-u", 1.65);
AcCritValues.put("60-4-0.05-l", 1.48);
AcCritValues.put("60-4-0.05-u", 1.69);
AcCritValues.put("60-5-0.05-l", 1.44);
AcCritValues.put("60-5-0.05-u", 1.73);
AcCritValues.put("60-6-0.05-l", 1.41);
AcCritValues.put("60-6-0.05-u", 1.77);
AcCritValues.put("60-2-0.01-l", 1.38);
AcCritValues.put("60-2-0.01-u", 1.45);
AcCritValues.put("60-3-0.01-l", 1.35);
AcCritValues.put("60-3-0.01-u", 1.48);
AcCritValues.put("60-4-0.01-l", 1.32);
AcCritValues.put("60-4-0.01-u", 1.52);
AcCritValues.put("60-5-0.01-l", 1.28);
AcCritValues.put("60-5-0.01-u", 1.56);
AcCritValues.put("60-6-0.01-l", 1.25);
AcCritValues.put("60-6-0.01-u", 1.60);

AcCritValues.put("80-2-0.05-l", 1.61);
AcCritValues.put("80-2-0.05-u", 1.66);
AcCritValues.put("80-3-0.05-l", 1.59);
AcCritValues.put("80-3-0.05-u", 1.69);
AcCritValues.put("80-4-0.05-l", 1.56);
AcCritValues.put("80-4-0.05-u", 1.72);
AcCritValues.put("80-5-0.05-l", 1.53);
AcCritValues.put("80-5-0.05-u", 1.74);
AcCritValues.put("80-6-0.05-l", 1.51);
AcCritValues.put("80-6-0.05-u", 1.77);
AcCritValues.put("80-2-0.01-l", 1.47);
AcCritValues.put("80-2-0.01-u", 1.52);
AcCritValues.put("80-3-0.01-l", 1.44);
AcCritValues.put("80-3-0.01-u", 1.54);
AcCritValues.put("80-4-0.01-l", 1.42);
AcCritValues.put("80-4-0.01-u", 1.57);
AcCritValues.put("80-5-0.01-l", 1.39);
AcCritValues.put("80-5-0.01-u", 1.60);
AcCritValues.put("80-6-0.01-l", 1.36);
AcCritValues.put("80-6-0.01-u", 1.62);

AcCritValues.put("100-2-0.05-l", 1.65);
AcCritValues.put("100-2-0.05-u", 1.69);
AcCritValues.put("100-3-0.05-l", 1.63);
AcCritValues.put("100-3-0.05-u", 1.72);
AcCritValues.put("100-4-0.05-l", 1.61);
AcCritValues.put("100-4-0.05-u", 1.74);
AcCritValues.put("100-5-0.05-l", 1.59);
AcCritValues.put("100-5-0.05-u", 1.76);
AcCritValues.put("100-6-0.05-l", 1.57);
AcCritValues.put("100-6-0.05-u", 1.78);
AcCritValues.put("100-2-0.01-l", 1.52);
AcCritValues.put("100-2-0.01-u", 1.56);
AcCritValues.put("100-3-0.01-l", 1.50);
AcCritValues.put("100-3-0.01-u", 1.58);
AcCritValues.put("100-4-0.01-l", 1.48);
AcCritValues.put("100-4-0.01-u", 1.60);
AcCritValues.put("100-5-0.01-l", 1.46);
AcCritValues.put("100-5-0.01-u", 1.63);
AcCritValues.put("100-6-0.01-l", 1.44);
AcCritValues.put("100-6-0.01-u", 1.65);
}

//endregion

//region White noise

/**  Get White Noise critical value 
@param n number of cases 
@param alpha significance level 
@return  critical value 
*/
//C# TO JAVA CONVERTER TODO TASK: The following line could not be converted:
public double GetWNCritValue(int n, double alpha)
{
return norm.qNorm(1.0 - alpha / 2.0, 0, 1) / Math.sqrt((double)n);
}

/**  White noise test (randomness) 
@param acfs autocorrelation function series 
@param alpha significance level 
@return  if the time series is white noise or not 
*/

public boolean IsWhiteNoise(List<Double> acfs, double alpha) {
int n = acfs.size();
double cv = GetWNCritValue(n, alpha);
double cvMax = GetWNCritValue(n, 0.01);
double ac;
double a = 0;
for (int i = 1;i < acfs.size();i++)
{
	ac = Math.abs(acfs.get(i));
	if (ac > cvMax)
	{
		return false;
	}
	if (ac > cv)
	{
		a++;
	}
}
if (a / (double) n > alpha)
{
	return false;
}
return true;
}

/**  Moving average test 
@param acfs autocorrelation function series 
@param alpha significance level 
@return  if the time series is a MA(q), returns q value. If not, returns -1 
*/
//C# TO JAVA CONVERTER TODO TASK: The following line could not be converted:
public int TestMovingAverage(List<Double> acfs, double alpha)
{
int n = acfs.size();
double cv = GetWNCritValue(n, alpha);
double ac;
for (int i = n - 1;i > 0;i--)
{
	ac = Math.abs(acfs.get(i));
	if (ac > cv)
	{
		return i + 1;
	}
	if (i == 0)
	{
		return i;
	}
}
return -1;
}

/**  Autorregresive test 
@param acfs autocorrelation function series 
@param alpha significance level 
@return  if the time series is a AR(p), returns p value. If not, returns -1 
*/

public int TestAutorregresive(List<Double> acfs, double alpha)
{
ArrayList<Double> pacfs = stat.PACF(acfs);
int n = pacfs.size();
double cv = GetWNCritValue(n, alpha);
double ac;
for (int i = n - 1;i > 0;i++)
{
	ac = Math.abs(pacfs.get(i));
	if (ac > cv)
	{
		return i + 1;
	}
	if (i == 0)
	{
		return i;
	}
}
return -1;
}

//endregion

//endregion

//region Variance stabilization

//C# TO JAVA CONVERTER TODO TASK: The following line could not be converted:
public double GetBoxCoxLambda(List<Double> serie, double lambdaInc, DotNetHelpers.RefObject<Double> maxLik)
{
maxLik.argValue = -Double.MAX_VALUE;
double maxLambda = -1;
double lik;
ArrayList<Double> transf;
double lambda = -3;
while (lambda <= 3)
{
	transf = (ArrayList<Double>) GetBoxCoxTransf(serie, lambda);
	lik = Likelihood(transf, lambda);
	if (lik > maxLik.argValue)
	{
		maxLik.argValue = lik;
		maxLambda = lambda;
	}
	lambda = Math.round((lambda + lambdaInc) * Math.pow(10, 2)) / Math.pow(10, 2);
}
return maxLambda;

}

//C# TO JAVA CONVERTER TODO TASK: The following line could not be converted:
public List<Double> GetBoxCoxTransf(List<Double> serie, double lambda)
{
ArrayList<Double> transf = new ArrayList<Double>();
for (double x : serie)
{
	transf.add(BoxCoxTransf(x, lambda));
}
return transf;
}

//C# TO JAVA CONVERTER TODO TASK: The following line could not be converted:
private double Likelihood(List<Double> serie, double lambda)
{
return this.ShapiroWilk(serie);
}

//C# TO JAVA CONVERTER TODO TASK: The following line could not be converted:
public double BoxCoxTransf(double x, double lambda)
{
if (lambda == 0)
{
	return Math.log(x);
}
else
{
	return (Math.pow(x, lambda) - 1) / lambda;
}
}

//C# TO JAVA CONVERTER TODO TASK: The following line could not be converted:
public double BoxCoxTransfInv(double y, double lambda)
{
if (lambda == 0)
{
	return Math.exp(y);
}
else
{
	return Math.pow(lambda * y + 1, 1 / lambda);
}
}

//endregion

}

