package ClassicalStat;


//region Imports

import java.util.*;
import DotNetHelpers.*;

//endregion


//DotNetHelpers_doc_comment_start  A new class for statistics calculation 
//DotNetHelpers_doc_comment_end
public class Statistics
{

    //region Constructor

    //DotNetHelpers_doc_comment_start  Constructor  
    //DotNetHelpers_doc_comment_end
    public Statistics()
    {
    }

    //endregion
    
    //region Public Methods
    //region Classic Statistics
    //region Basic Functions

    //DotNetHelpers_doc_comment_start  Mean of the collection 
    //DotNetHelpers_doc_comment_body @param values collection of values 
    //DotNetHelpers_doc_comment_body @return  returns the mean 
    //DotNetHelpers_doc_comment_end
    public final double Mean(List<Double> values)
    {
        if(values.isEmpty())
        {
            return 0.0;
        }

        double tot = 0.0;
        for(int i=0;i<values.size();i++)
        {
            tot += values.get(i);
        }
        return tot/values.size();
    }

    //DotNetHelpers_doc_comment_start  Variance of the collection 
    //DotNetHelpers_doc_comment_body @param values collection of values 
    //DotNetHelpers_doc_comment_body @return  returns the variance 
    //DotNetHelpers_doc_comment_end
    public final double Variance(List<Double> values)
    {
        int n = values.size();
        if(n == 1)
        {
            return 0;
        }
        double Sum = 0.0, SumSquares = 0.0;

        for(int i=0;i<n;i++)
        {
            Sum += values.get(i);
            SumSquares += values.get(i) * values.get(i);
        }
        return (n * SumSquares - (Sum * Sum)) / (n*n - 1);
    }

    //DotNetHelpers_doc_comment_start  Range of a collection 
    //DotNetHelpers_doc_comment_body @param values collection of values 
    //DotNetHelpers_doc_comment_body @return  returns the range 
    //DotNetHelpers_doc_comment_end
    public final double Range(List<Double> values)
    {
        double min = Double.MAX_VALUE;
        double max = -Double.MAX_VALUE;
        for (double val : values)
        {
            if (val < min)
            {
                min = val;
            }
            if (val > max)
            {
                max = val;
            }
        }
        return max - min;
    }

    //DotNetHelpers_doc_comment_start  Standard deviation of the collection 
    //DotNetHelpers_doc_comment_body @param values collection of values 
    //DotNetHelpers_doc_comment_body @return  returns the standard deviation 
    //DotNetHelpers_doc_comment_end
    public final double StDev(List<Double> values)
    {
        return Math.sqrt(Variance(values));
    }

    //DotNetHelpers_doc_comment_start  Variation coefficent of the collection 
    //DotNetHelpers_doc_comment_body @param values collection of values 
    //DotNetHelpers_doc_comment_body @return  returns the variation coefficent 
    //DotNetHelpers_doc_comment_end
    public final double VarCoeff(List<Double> values)
    {
        if(values.size() <= 1)
        {
            return 0;
        }
        double mean = Mean(values);
        if(mean == 0)
        {
            return 1;
        }
        double varCoeff = Math.abs(StDev(values)/ mean);
        return varCoeff;
    }

    //DotNetHelpers_doc_comment_start  Variance from pre-calculated sums 
    //DotNetHelpers_doc_comment_body @param count quantity of elements 
    //DotNetHelpers_doc_comment_body @param sum sum of elements 
    //DotNetHelpers_doc_comment_body @param sqSum sum of squares 
    //DotNetHelpers_doc_comment_body @return  the variance 
    //DotNetHelpers_doc_comment_end
    public final double Variance(int count, double sum, double sqSum)
    {
        if(count <= 1 || sum == 0)
        {
            return 0.0;
        }
        return (count * sqSum - (sum * sum)) / (count*count - 1);
    }


    //DotNetHelpers_doc_comment_start  Standard deviation from pre-calculated sums 
    //DotNetHelpers_doc_comment_body @param count quantity of elements 
    //DotNetHelpers_doc_comment_body @param sum sum of elements 
    //DotNetHelpers_doc_comment_body @param sqSum sum of squares 
    //DotNetHelpers_doc_comment_body @return  the standard deviation 
    //DotNetHelpers_doc_comment_end
    public final double StDev(int count, double sum, double sqSum)
    {
        double variance = Variance(count, sum, sqSum);
        if(variance < 0.01)
        {
            return 0.0;
        }
        return Math.sqrt(variance);
    }

    //DotNetHelpers_doc_comment_start  Standard error of the mean 
    //DotNetHelpers_doc_comment_body @param values list of values 
    //DotNetHelpers_doc_comment_body @return  the standard error 
    //DotNetHelpers_doc_comment_end
    public final double StdError(List<Double> values)
    {
        if(values.isEmpty())
        {
            throw new RuntimeException("Strings.Empty_list");
        }
        return StDev(values)/Math.sqrt(values.size());
    }

    //DotNetHelpers_doc_comment_start  Mean Square Error between two collections 
    //DotNetHelpers_doc_comment_body @param data1 first collection 
    //DotNetHelpers_doc_comment_body @param data2 second collection 
    //DotNetHelpers_doc_comment_body @return  the MSE value 
    //DotNetHelpers_doc_comment_end
    public final double MSError(List<Double> data1, List<Double> data2)
    {
        if(data1.size() != data2.size())
        {
            throw new RuntimeException("Strings.s1_must_have_same_size_as_s2");
        }

        double tot = 0.0;
        for(int i=0;i<data1.size();i++)
        {
            tot += Math.pow((data1.get(i) - data2.get(i)),2);
        }
        return Math.sqrt(tot)/(double)data1.size();
    }

    //DotNetHelpers_doc_comment_start  Mean Square Error from two collections with a percentage for checking 
    //DotNetHelpers_doc_comment_body @param data1 the first collection 
    //DotNetHelpers_doc_comment_body @param data2 the second collection 
    //DotNetHelpers_doc_comment_body @param percCheck percentage for checking 
    //DotNetHelpers_doc_comment_body @return  the mse value 
    //DotNetHelpers_doc_comment_end
    public final double MSError(List<Double> data1, List<Double> data2, double percCheck)
    {
        if(data1.size() > data2.size())
        {
            List<Double> aux = data1;
            data1 = data2;
            data2 = aux;
        }
        int firstIndex = data2.size() - (int)(percCheck * 0.01 * data2.size());
        double tot = 0.0;
        for(int i=firstIndex;i<data1.size();i++)
        {
            tot += Math.pow((data1.get(i) - data2.get(i)), 2);
        }
        for(int i=data1.size();i<data2.size();i++)
        {
            tot += Math.pow(data2.get(i), 2);
        }
        return Math.sqrt(tot)/(double)(data2.size() - firstIndex);
    }

    //DotNetHelpers_doc_comment_start  Covariace 
    //DotNetHelpers_doc_comment_body @param data1 the first collection 
    //DotNetHelpers_doc_comment_body @param data2 the second collection 
    //DotNetHelpers_doc_comment_body @return  covariance value 
    //DotNetHelpers_doc_comment_end
    public final double Cov(List<Double> data1, List<Double> data2)
    {
        if(data1.size() != data2.size())
        {
            throw new RuntimeException("Strings.Sets_must_have_the_same_size");
        }

        double mean1 = Mean(data1);
        double mean2 = Mean(data2);
        double sumProds = 0;
        for(int i = 0;i < data1.size();i++)
        {
            sumProds += ((data1.get(i) - mean1) * (data2.get(i) - mean2));
        }
        return sumProds/data1.size();
    }

    //DotNetHelpers_doc_comment_start  
    //DotNetHelpers_doc_comment_body Pearson R correlation coefficient
    //DotNetHelpers_doc_comment_body Range [-1,1] : 1 full direct correlation, -1 full inverse correlation, 0 no correlation.
    //DotNetHelpers_doc_comment_body 
    //DotNetHelpers_doc_comment_body @param data1 the first collection 
    //DotNetHelpers_doc_comment_body @param data2 the second collection 
    //DotNetHelpers_doc_comment_body @return  r value 
    //DotNetHelpers_doc_comment_end
    public final double R(List<Double> data1, List<Double> data2)
    {
        if(StDev(data1) * StDev(data2) == 0)
        {
            return 0;
        }
        return Cov(data1, data2) / (StDev(data1) * StDev(data2));
    }

    //DotNetHelpers_doc_comment_start  R squared determination coefficient 
    //DotNetHelpers_doc_comment_body @param data1 the first collection 
    //DotNetHelpers_doc_comment_body @param data2 the second collection 
    //DotNetHelpers_doc_comment_body @return  r2 value 
    //DotNetHelpers_doc_comment_end
    public final double R2(List<Double> data1, List<Double> data2)
    {
        double r2 = Math.pow(R(data1, data2), 2);
        if (r2 < 0.0001)
        {
            r2 = 0;
        }
        return r2;
    }

    //DotNetHelpers_doc_comment_start  Adjusted R squared for non-nested model comparison 
    //DotNetHelpers_doc_comment_body @param data1 the first collection 
    //DotNetHelpers_doc_comment_body @param data2 the second collection 
    //DotNetHelpers_doc_comment_body @param n size of collection 
    //DotNetHelpers_doc_comment_body @param p number of parameters 
    //DotNetHelpers_doc_comment_body @return  adjusted r2 value 
    //DotNetHelpers_doc_comment_end
    public final double AdjR2(List<Double> data1, List<Double> data2, double p)
    {
        if(data1.size() != data2.size())
        {
            throw new RuntimeException("Strings.data1_and_data2_must_have_the_same_size");
        }
        int n = data1.size();
        double r2 = R2(data1, data2);
        return 1 - (1-r2) * (n-1)/(n-p-1);
    }

    //DotNetHelpers_doc_comment_start  sum of squares 
    //DotNetHelpers_doc_comment_body @param data data set 
    //DotNetHelpers_doc_comment_body @return  ss 
    //DotNetHelpers_doc_comment_end
    public final double SumSq(List<Double> data)
    {
        double sse = 0.0;
        for (int i = 0; i < data.size(); i++)
        {
            sse = sse + data.get(i) * data.get(i);
        }
        return sse;
    }

    //endregion

    //region Interval functions

    public final double Mean(List<Double> values, int ini, int end)
    {
        if (values.isEmpty())
        {
            return 0.0;
        }
        if (ini < 0 || ini >= end || end > values.size())
        {
            throw new RuntimeException("Error");
        }

        double tot = 0.0;
        for (int i = ini; i <= end; i++)
        {
            tot += values.get(i);
        }
        return tot / (end - ini + 1);
    }

    public final double Variance(List<Double> values, int ini, int end)
    {
        if (ini < 0 || ini >= end || end > values.size())
        {
            throw new RuntimeException("Error");
        }

        int n = end - ini + 1;
        if (n == 1)
        {
            return 0;
        }
        double Sum = 0.0, SumSquares = 0.0;

        for (int i = ini; i <= end; i++)
        {
            Sum += values.get(i);
            SumSquares += values.get(i) * values.get(i);
        }
        double var = (SumSquares - (Sum * Sum) / n) / (n - 1);
        if (var < 0)
        {
            var = 0;
        }
        return var;
    }

    public final double Range(List<Double> values, int ini, int end)
    {
        if (ini < 0 || ini >= end || end > values.size())
        {
            throw new RuntimeException("Error");
        }

        double min = Double.MAX_VALUE;
        double max = -Double.MAX_VALUE;
        for(int i=ini;i<=end;i++)
        {
            if (values.get(i) < min)
            {
                min = values.get(i);
            }
            if (values.get(i) > max)
            {
                max = values.get(i);
            }
        }
        return max - min;
    }

    public final double StDev(List<Double> values, int ini, int end)
    {
        return Math.sqrt(Variance(values, ini, end));
    }

    public final double Cov(List<Double> data1, List<Double> data2, int ini, int end)
    {
        if (data1.size() != data2.size())
        {
            throw new RuntimeException("Strings.Sets_must_have_the_same_size");
        }
        if (ini < 0 || ini >= end || end > data1.size())
        {
            throw new RuntimeException("Error");
        }

        double mean1 = Mean(data1, ini, end);
        double mean2 = Mean(data2, ini, end);
        double sumProds = 0;
        for (int i = ini; i <= end; i++)
        {
            sumProds += ((data1.get(i) - mean1) * (data2.get(i) - mean2));
        }
        return sumProds / (end - ini + 1);
    }

    public final double R(List<Double> X, List<Double> Y, int ini, int end)
    {
        double m1 = Mean(X, ini, end);
        double m2 = Mean(Y, ini, end);
        double s1 = StDev(X, ini, end);
        double s2 = StDev(Y, ini, end);
        if (s1 * s2 == 0)
        {
            if(m1 == m2)
            {
                return 1;
            }
            else
            {
                return 0;
            }
        }

        double n = end - ini + 1;
        double sX = 0;
        double sY = 0;
        double sX2 = 0;
        double sY2 = 0;
        double sXY = 0;
        for (int i = ini; i <= end; i++)
        {
            sX += X.get(i);
            sX2 += (X.get(i) * X.get(i));
            sY += Y.get(i);
            sY2 += (Y.get(i) * Y.get(i));
            sXY += (X.get(i) * Y.get(i));

        }

        double r = (n * sXY - sX * sY) / Math.sqrt((n * sX2 - Math.pow(sX, 2)) * (n * sY2 - Math.pow(sY, 2)));
        if (r < 0)
        {
            r = 0;
        }
        return r;
    }

    public final double R2(List<Double> data1, List<Double> data2, int ini, int end)
    {
        return Math.pow(R(data1, data2, ini, end), 2);
    }

    public final ArrayList<Double> R2Mw(List<Double> data1, List<Double> data2, int mw)
    {
        if (data1.size() != data2.size())
        {
            throw new RuntimeException("Strings.Sets_must_have_the_same_size");
        }
        ArrayList<Double> r2List = new ArrayList<Double>();
        double r2;
        for (int i = 0; i < data1.size()-mw; i++)
        {
            r2 = R2(data1, data2, i, i + mw);
            r2List.add(r2);
        }
        return r2List;
    }

    public final ArrayList<Double> R2ToEnd(List<Double> data1, List<Double> data2)
    {
        if (data1.size() != data2.size())
        {
            throw new RuntimeException("Strings.Sets_must_have_the_same_size");
        }
        ArrayList<Double> r2List = new ArrayList<Double>();
        double r2;
        for (int i = 0; i < data1.size()-1; i++)
        {
            r2 = R2(data1, data2, i, data1.size());
            r2List.add(r2);
        }
        return r2List;
    }

    public final double R(List<Double> data1, List<Double> data2, int ini, int end, int lag)
    {
        if(ini + lag >= end)
        {
            return 0;
        }
        ArrayList<Double> data1Lag = new ArrayList<Double>(data1);
        DotNetHelpers.DoubleLists.RemoveRange(data1Lag, 0, lag + 0);
        ArrayList<Double> data2Lag = new ArrayList<Double>(data2);
        DotNetHelpers.DoubleLists.RemoveRange(data2Lag, data2.size()-lag, lag + data2.size()-lag);
        return R(data1Lag, data2Lag, ini, end-lag);
    }

    public final double R2(List<Double> data1, List<Double> data2, int ini, int end, int lag)
    {
        return Math.pow(R(data1, data2, ini, end, lag), 2);
    }

    public final ArrayList<Double> R2Mw(List<Double> data1, List<Double> data2, int mw, int lag)
    {
        if (data1.size() != data2.size())
        {
            throw new RuntimeException("Strings.Sets_must_have_the_same_size");
        }
        ArrayList<Double> data1Lag = new ArrayList<Double>(data1);
        DotNetHelpers.DoubleLists.RemoveRange(data1Lag, lag, lag+0);
        ArrayList<Double> data2Lag = new ArrayList<Double>(data2);
        DotNetHelpers.DoubleLists.RemoveRange(data2Lag,data2.size() - lag, lag + data2.size() - lag);
        ArrayList<Double> r2List = new ArrayList<Double>();
        double r2;
        for (int i = 0; i < data1Lag.size() - mw; i++)
        {
            r2 = R2(data1Lag, data2Lag, i, i + mw);
            r2List.add(r2);
        }
        return r2List;
    }

    public final int FirstR2AboveThres(ArrayList<Double> r2List, double threshold)
    {
        if(r2List.get(r2List.size()-1).compareTo(threshold) < 0)
        {
            return -1;
        }
        for (int i = r2List.size() - 2; i > 0; i--)
        {
            if (r2List.get(i).compareTo(threshold) < 0)
            {
                return i + 1;
            }
        }
        return 0;
    }

    //endregion

    //region Memory functions

    public final void CalcSums(List<Double> serie, List<Double> weights, DotNetHelpers.RefObject<Double> count, DotNetHelpers.RefObject<Double> sum, DotNetHelpers.RefObject<Double> sumSq)
    {
        count.argValue = 0.0;
        sum.argValue = 0.0;
        sumSq.argValue = 0.0;
        for (int i = 0; i < serie.size(); i++)
        {
            sum.argValue = sum.argValue + serie.get(i) * weights.get(i);
            sumSq.argValue = sumSq.argValue + serie.get(i) * serie.get(i) * weights.get(i);
            count.argValue = count.argValue + weights.get(i);
        }
    }

    public final double Mean(double count, double sum)
    {
        return sum / count;
    }

    public final double Variance(double count, double sum, double sumSq)
    {
        return (count * sumSq - (sum * sum)) / (count * count - 1);
    }

    public final double StDev(double count, double sum, double sumSq)
    {
        double variance = Variance(count, sum, sumSq);
        if (variance < 0.01)
        {
            return 0.0;
        }
        return Math.sqrt(variance);
    }

    public final double WCov(List<Double> data1, List<Double> data2, List<Double> weights)
    {
        if (data1.size() != data2.size())
        {
            throw new RuntimeException("Strings.Sets_must_have_the_same_size");
        }
        double mean1 = Mean(data1);
        double mean2 = Mean(data2);
        double sumProds = 0;
        double count = 0;
        for (int i = 0; i < data1.size(); i++)
        {
            sumProds = sumProds + (data1.get(i) - mean1) * (data2.get(i) - mean2) * weights.get(i);
            count = count + weights.get(i);
        }
        return sumProds / count;
    }


    public final double WR(double sd1, double sd2, double cov)
    {
        if (sd1 * sd2 == 0)
        {
            return 0;
        }
        return cov / (sd1 * sd2);
    }

    public final double WR2(double sd1, double sd2, double cov)
    {
        double r2 = Math.pow(WR(sd1, sd2, cov), 2);
        if (r2 < 0.0001)
        {
            r2 = 0;
        }
        return r2;
    }

    public final double WR2(List<Double> data1, List<Double> data2, List<Double> weights)
    {
        double count1 = 0;
        double sum1 = 0;
        double sumSq1 = 0;
        double count2 = 0;
        double sum2 = 0;
        double sumSq2 = 0;
        DotNetHelpers.RefObject<Double> tempRef_count1 = new DotNetHelpers.RefObject<Double>(count1);
        DotNetHelpers.RefObject<Double> tempRef_sum1 = new DotNetHelpers.RefObject<Double>(sum1);
        DotNetHelpers.RefObject<Double> tempRef_sumSq1 = new DotNetHelpers.RefObject<Double>(sumSq1);
        CalcSums(data1, weights, tempRef_count1, tempRef_sum1, tempRef_sumSq1);
    sumSq1 = tempRef_sumSq1.argValue;
    sum1 = tempRef_sum1.argValue;
    count1 = tempRef_count1.argValue;
        DotNetHelpers.RefObject<Double> tempRef_count2 = new DotNetHelpers.RefObject<Double>(count2);
        DotNetHelpers.RefObject<Double> tempRef_sum2 = new DotNetHelpers.RefObject<Double>(sum2);
        DotNetHelpers.RefObject<Double> tempRef_sumSq2 = new DotNetHelpers.RefObject<Double>(sumSq2);
        CalcSums(data2, weights, tempRef_count2, tempRef_sum2, tempRef_sumSq2);
    sumSq2 = tempRef_sumSq2.argValue;
    sum2 = tempRef_sum2.argValue;
    count2 = tempRef_count2.argValue;
        double sd1 = StDev(count1, sum1, sumSq1);
        double sd2 = StDev(count2, sum2, sumSq2);
        double cov = WCov(data1, data2, weights);
        return WR2(sd1, sd2, cov);
    }

    public final ArrayList<Double> CalculateR2(List<Double> serie1, List<Double> serie2, int movingWindow)
    {
        ArrayList<Double> r2s = new ArrayList<Double>();
        int minCount = Math.min(serie1.size(), serie2.size());
        for (int i = 0; i < minCount - movingWindow - 1; i++)
        {
            double r2 = R2(serie1, serie2, i, (i + movingWindow - 1));
            r2s.add(r2);
        }
        return r2s;
    }

    public final double WeightedR2(List<Double> r2s, double[] weights)
    {
        double sum = 0;
        for (int i = 0; i < r2s.size(); i++)
        {
            sum = sum + r2s.get(i) * weights[i];
        }
        return sum / r2s.size();
    }

    //endregion

    //region Basic functions in arrays

    //DotNetHelpers_doc_comment_start  Calculate a mean in a list ignoring the datum of the index 
    //DotNetHelpers_doc_comment_body @param datos IList of data 
    //DotNetHelpers_doc_comment_body @param indice the index to ignore 
    //DotNetHelpers_doc_comment_body @return  the calculated mean 
    //DotNetHelpers_doc_comment_end
    public final double MeanWithout(List<Double> datos, int indice)
    {
        double total = 0.0;
        for(int i=0;i<datos.size();i++)
        {
            if(i != indice)
            {
                total += datos.get(i);
            }
        }
        return total / (datos.size() - 1);
    }

    //DotNetHelpers_doc_comment_start  Calculate a mean in a list ignoring the datum of the index 
    //DotNetHelpers_doc_comment_body @param datos IList of data 
    //DotNetHelpers_doc_comment_body @param indice the index to ignore 
    //DotNetHelpers_doc_comment_body @return  the calculated mean 
    //DotNetHelpers_doc_comment_end
    public final double StDevWithout(List<Double> datos, int indice)
    {
        double total = 0.0;
        double tot2 = 0.0;
        for(int i=0;i<datos.size();i++)
        {
            if(i != indice)
            {
                total += datos.get(i);
                tot2 += Math.pow(datos.get(i), 2);
            }
        }
        return Math.sqrt(((datos.size() - 1) * tot2 - Math.pow(total, 2)) / ((datos.size() - 1) * (datos.size() - 2)));

    }

    //DotNetHelpers_doc_comment_start  Calculate a bias in a list ignoring the datum of the index 
    //DotNetHelpers_doc_comment_body @param datos IList of data 
    //DotNetHelpers_doc_comment_body @param indice the index to ignore 
    //DotNetHelpers_doc_comment_body @return  the calculated bias 
    //DotNetHelpers_doc_comment_end
    public final double BiasWithout(List<Double> datos, int indice)
    {
        double mean = MeanWithout(datos, indice);
        return Math.abs(datos.get(indice) - mean);
    }

    //DotNetHelpers_doc_comment_start  Calculate a mean of an interval within a list 
    //DotNetHelpers_doc_comment_body @param datos data for calculation 
    //DotNetHelpers_doc_comment_body @param indice index of initial period 
    //DotNetHelpers_doc_comment_body @param r radium 
    //DotNetHelpers_doc_comment_body @return  the calculated mean 
    //DotNetHelpers_doc_comment_end
    public final double IntervalMean(List<Double> datos, int indice, int r)
    {
        double total = 0.0;
        int n = 0;
        for(int i=indice-r;i<=indice+r;i++)
        {
            if(i != indice && i >= 0 && i < datos.size())
            {
                total += datos.get(i);
                n++;
            }
        }
        return total / n;
    }

    //DotNetHelpers_doc_comment_start  Calculate a bias of an interval within a list 
    //DotNetHelpers_doc_comment_body @param datos data for calculation 
    //DotNetHelpers_doc_comment_body @param indice index of initial period 
    //DotNetHelpers_doc_comment_body @param r radium 
    //DotNetHelpers_doc_comment_body @return  the calculated bias 
    //DotNetHelpers_doc_comment_end
    public final double IntervalBias(List<Double> datos, int indice, int r)
    {
        double mean = IntervalMean(datos, indice, r);
        return Math.abs(datos.get(indice) - mean);
    }
    //endregion
    //endregion
    //region Robust statistics
    //region Trimmed statistics 

    //DotNetHelpers_doc_comment_start  Trim a list of values 
    //DotNetHelpers_doc_comment_body @param values list of values 
    //DotNetHelpers_doc_comment_body @param cutoff cut off in each boundaries 
    //DotNetHelpers_doc_comment_body @param sort if values need to be sorted 
    //DotNetHelpers_doc_comment_body @return  trimmed list of values 
    //DotNetHelpers_doc_comment_end
    public final ArrayList<Double> Trim(ArrayList<Double> values, double cutoff, boolean sort)
    {
        if(sort)
        {
            Collections.sort(values);
        }
        int half = (int)((double)values.size()/2);
        int nCutoff = (int)java.lang.Math.round((double)half * cutoff);
        int nValues = values.size() - nCutoff*2;
        if(nValues == 0)
        {
            return values;
        }
        ArrayList<Double> trimValues = new ArrayList<Double>();
        for(int i = nCutoff;i < nValues-nCutoff;i++)
        {
            trimValues.add(values.get(i));
        }
        return trimValues;
    }

    //DotNetHelpers_doc_comment_start  Trimmed Range 
    //DotNetHelpers_doc_comment_body @param values list of values 
    //DotNetHelpers_doc_comment_body @param cutoff cut off in each boundaries 
    //DotNetHelpers_doc_comment_body @param sort if values need to be sorted 
    //DotNetHelpers_doc_comment_body @return  trimmed range 
    //DotNetHelpers_doc_comment_end
    public final double TrimRange(ArrayList<Double> values, double cutoff, boolean sort)
    {
        ArrayList<Double> trimValues = Trim(values, cutoff, sort);
        if (values.isEmpty())
        {
            return 0.0;
        }
        return Range(trimValues);
    }

    public final double TrimRange(ArrayList<Double> values, int ini, int end, double cutoff, boolean sort)
    {
        ArrayList<Double> trimValues = Trim(values, cutoff, sort);
        if (values.isEmpty())
        {
            return 0.0;
        }
        return Range(trimValues);
    }


    //DotNetHelpers_doc_comment_start  Trimmed Mean 
    //DotNetHelpers_doc_comment_body @param values list of values 
    //DotNetHelpers_doc_comment_body @param cutoff cut off in each boundaries 
    //DotNetHelpers_doc_comment_body @param sort if values need to be sorted 
    //DotNetHelpers_doc_comment_body @return  trimmed mean 
    //DotNetHelpers_doc_comment_end
    public final double TrimMean(ArrayList<Double> values, double cutoff, boolean sort)
    {
        ArrayList<Double> trimValues = Trim(values, cutoff, sort);
        if(values.isEmpty())
        {
            return 0.0;
        }
        return Mean(trimValues);
    }

    //DotNetHelpers_doc_comment_start  Trimmed Variance 
    //DotNetHelpers_doc_comment_body @param values list of values 
    //DotNetHelpers_doc_comment_body @param cutoff cut off in each boundaries 
    //DotNetHelpers_doc_comment_body @param sort if values need to be sorted 
    //DotNetHelpers_doc_comment_body @return  trimmed variance 
    //DotNetHelpers_doc_comment_end
    public final double TrimVar(ArrayList<Double> values, double cutoff, boolean sort)
    {
        ArrayList<Double> trimValues = Trim(values, cutoff, sort);
        if(trimValues.size() <= 1)
        {
            return 0.0;
        }
        return Variance(trimValues);
    }

    //DotNetHelpers_doc_comment_start  Trimmed Standard deviation 
    //DotNetHelpers_doc_comment_body @param values list of values 
    //DotNetHelpers_doc_comment_body @param cutoff cut off in each boundaries 
    //DotNetHelpers_doc_comment_body @param sort if values need to be sorted 
    //DotNetHelpers_doc_comment_body @return  trimmed standard deviation 
    //DotNetHelpers_doc_comment_end
    public final double TrimSd(ArrayList<Double> values, double cutoff, boolean sort)
    {
        ArrayList<Double> trimValues = Trim(values, cutoff, sort);
        if(trimValues.size() <= 1)
        {
            return 0.0;
        }
        return StDev(trimValues);
    }

    //DotNetHelpers_doc_comment_start  Trimmed Mean 
    //DotNetHelpers_doc_comment_body @param values list of values 
    //DotNetHelpers_doc_comment_body @param rec cut off in each boundaries 
    //DotNetHelpers_doc_comment_body @return  trimmed mean 
    //DotNetHelpers_doc_comment_end
    public final double TrimMeanNonZero(ArrayList<Double> values, double rec)
    {
        int mitad = (int)((double)values.size()/2);
        int nRec = (int)java.lang.Math.round((double)mitad * rec);
        int nValues = 0;
        double sum = 0.0;
        for(int i = nRec;i < values.size()-nRec;i++)
        {
            if( ! values.get(i).equals(0))
            {
                sum += values.get(i);
                nValues++;
            }
        }
        if(nValues == 0)
        {
            return 0;
        }
        return sum/nValues;
    }

    //endregion

    //region Windsored statistics 

    //DotNetHelpers_doc_comment_start  Winsored Mean 
    //DotNetHelpers_doc_comment_body @param values list of values 
    //DotNetHelpers_doc_comment_body @param rec cut off in each boundaries 
    //DotNetHelpers_doc_comment_body @return  winsored mean 
    //DotNetHelpers_doc_comment_end
    public final double WinMean(double[] values, double rec)
    {
        int mitad = (int)((double)values.length/2);
        int nRec = (int)java.lang.Math.round((double)mitad * (1 -rec));
        int nValues = values.length;
        if(nValues == 0)
        {
            return 0.0;
        }
        double sum = 0.0;
        for(int i = 0;i<nRec;i++)
        {
            sum += values[nRec];
        }
        for(int i = nRec;i<nValues-nRec;i++)
        {
            sum += values[i];
        }
        for(int i= nValues-nRec;i<nValues;i++)
        {
            sum += values[nValues-nRec-1];
        }
        return sum/nValues;
    }

    //DotNetHelpers_doc_comment_start  Windsor a list of values 
    //DotNetHelpers_doc_comment_body @param values list of values 
    //DotNetHelpers_doc_comment_body @param cutoff cut off in each boundaries 
    //DotNetHelpers_doc_comment_body @param sort if values need to be sorted 
    //DotNetHelpers_doc_comment_body @return  windsored list of values 
    //DotNetHelpers_doc_comment_end
    public final ArrayList<Double> Windsor(ArrayList<Double> values, double cutoff, boolean sort)
    {
        ArrayList<Double> valSort = values;
        if(sort)
        {
            valSort = new ArrayList<Double>(values);
            Collections.sort(valSort);
        }
        int half = (int)((double)valSort.size()/2);
        int nCutoff = (int)java.lang.Math.round((double)half * cutoff);
        int nValues = valSort.size() - nCutoff*2;
        if(nValues == 0)
        {
            return values;
        }
        double lwr = valSort.get(nCutoff);
        double upr = valSort.get(valSort.size()-nCutoff-1);

        ArrayList<Double> windValues = new ArrayList<Double>();
        for(int i=0;i<values.size();i++)
        {
            if (values.get(i).compareTo(lwr) < 0)
            {
                windValues.add(lwr);
            }
            else if (values.get(i).compareTo(upr) > 0)
            {
                windValues.add(upr);
            }
            else
            {
                windValues.add(values.get(i));
            }
        }
        return windValues;
    }

    //DotNetHelpers_doc_comment_start  Windsored Mean 
    //DotNetHelpers_doc_comment_body @param values list of values 
    //DotNetHelpers_doc_comment_body @param cutoff cut off in each boundaries 
    //DotNetHelpers_doc_comment_body @param sort if values need to be sorted 
    //DotNetHelpers_doc_comment_body @return  windsored mean 
    //DotNetHelpers_doc_comment_end
    public final double WindMean(ArrayList<Double> values, double cutoff, boolean sort)
    {
        ArrayList<Double> windValues = Windsor(values, cutoff, sort);
        if(values.isEmpty())
        {
            return 0.0;
        }
        return Mean(windValues);
    }

    public final ArrayList<Double> GetList(ArrayList<Double> values, int ini, int end)
    {
        ArrayList<Double> periodList = new ArrayList<Double>();
        for (int i = ini; i < end; i++)
        {
            periodList.add(values.get(i));
        }
        return periodList;
    }

    //DotNetHelpers_doc_comment_start  Windsored variance 
    //DotNetHelpers_doc_comment_body @param values list of values 
    //DotNetHelpers_doc_comment_body @param cutoff cut off in each boundaries 
    //DotNetHelpers_doc_comment_body @param sort if values need to be sorted 
    //DotNetHelpers_doc_comment_body @return  windsored variance 
    //DotNetHelpers_doc_comment_end
    public final double WindVar(ArrayList<Double> values, double cutoff, boolean sort)
    {
        ArrayList<Double> windValues = Windsor(values, cutoff, sort);
        if(values.size() <= 1)
        {
            return 0.0;
        }
        return Variance(windValues);
    }

    //DotNetHelpers_doc_comment_start  Windsored standard deviation 
    //DotNetHelpers_doc_comment_body @param values list of values 
    //DotNetHelpers_doc_comment_body @param cutoff cut off in each boundaries 
    //DotNetHelpers_doc_comment_body @param sort if values need to be sorted 
    //DotNetHelpers_doc_comment_body @return  windsored standard deviation 
    //DotNetHelpers_doc_comment_end
    public final double WindSd(ArrayList<Double> values, double cutoff, boolean sort)
    {
        ArrayList<Double> windValues = Windsor(values, cutoff, sort);
        if(values.size() <= 1)
        {
            return 0.0;
        }
        return StDev(windValues);
    }

    //endregion

    //endregion

    //region Entropy

    //DotNetHelpers_doc_comment_start  Calculate entropy of a collection of values 
    //DotNetHelpers_doc_comment_body @param values the collection of values 
    //DotNetHelpers_doc_comment_body @return  the entrophy value 
    //DotNetHelpers_doc_comment_end
    public final double CalcEntropy(List<Double> values)
    {
        NormalDist nd = new NormalDist();
        double mean = Mean(values);
        double stDev = StDev(values);
        ArrayList<Double> stdValues = new ArrayList<Double>();
        if(stDev == 0)
        {
            return 0.0;
        }
        for(double val : values)
        {
            stdValues.add((val - mean)/stDev);
        }
        double p, info;
        double entropy = 0.0;
        for(double val : stdValues)
        {
            p = nd.pNorm(val - 0.5, val + 0.5, 0, 1);
            info = Math.log(1.0/p);
            entropy += p * info;
        }
        return entropy;
    }
    //endregion
    //region Distributions
    //region t student

    //DotNetHelpers_doc_comment_start 
    //DotNetHelpers_doc_comment_body student t distribution
    //DotNetHelpers_doc_comment_body 
    //DotNetHelpers_doc_comment_body @param t t value 
    //DotNetHelpers_doc_comment_body @param df degrees of freedom 
    //DotNetHelpers_doc_comment_body @return  probability 
    //DotNetHelpers_doc_comment_end
    public final double pt(double t, int df)
    {
        if(df <=0)
        {
            throw new RuntimeException("Strings.Degrees_Of_Freedom_Must_Be_0");
        }
        double epsilon = 0.00001;
        double result = 0;
        double x = 0;
        double rk = 0;
        double z = 0;
        double f = 0;
        double tz = 0;
        double p = 0;
        double xsqk = 0;
        int j = 0;

        if(t == 0.0)
        {
            return 0.5;
        }
        if(t < -2.0)
        {
            rk = df;
            z = rk/(rk + t*t);
            result = 0.5 * Stat.betai(0.5 * rk, 0.5, z);
            return result;
        }
        if(t < 0)
        {
            x = -t;
        }
        else
        {
            x = t;
        }
        rk = df;
        z = 1.0+x*x/rk;
        if(df%2!=0)
        {
            xsqk = x/Math.sqrt(rk);
            p = Math.atan(xsqk);
            if(df>1)
            {
                f = 1.0;
                tz = 1.0;
                j = 3;
                while(j<=df-2 & (double)(tz/f) > epsilon)
                {
                    tz = tz*((j-1)/(z*j));
                    f = f+tz;
                    j = j+2;
                }
                p = p + f*xsqk/z;
            }
            p = p*2.0/Math.PI;
        }
        else
        {
            f = 1.0;
            tz = 1.0;
            j = 2;
            while(j<=df-2 & (double)(tz/f) > epsilon)
            {
                tz = tz*((j-1)/(z*j));
                f = f+tz;
                j = j+2;
            }
            p = f*x / Math.sqrt(z*rk);
        }
        if(t < 0)
        {
            p = -p;
        }
        result = 0.5 + 0.5*p;
        return result;
    }

    //DotNetHelpers_doc_comment_start  Inverse student t distribution 
    //DotNetHelpers_doc_comment_body @param p probability 
    //DotNetHelpers_doc_comment_body @param k degrees of freedom 
    //DotNetHelpers_doc_comment_body @return  quantil 
    //DotNetHelpers_doc_comment_end
    public final double qt(double p, int k)
    {
        if(p <= 0 || p > 1)
        {
            throw new RuntimeException("Strings.Probability_must_be_between_0_and_1");
        }
        double maxReal = 1000000;
        double result = 0;
        double t = 0;
        double rk = 0;
        double z = 0;
        int rflg = 0;

        rk = k;
        if(p > 0.25 && p < 0.75)
        {
            if(p == 0.5)
            {
                result = 0;
                return result;
            }
            z = 1.0-2.0 * p;
            z = Stat.betacf(0.5, 0.5*rk, Math.abs(z));
            t = Math.sqrt(rk*z/(1.0-z));
            if(p < 0.5)
            {
                t = -t;
            }
            result = t;
            return result;
        }
        rflg = -1;
        if(p >= 0.5)
        {
            p = 1.0-p;
            rflg = 1;
        }
        z = Stat.betacf(0.5*rk, 0.5, 2.0*p);
        if((double)(maxReal*z)< rk)
        {
            result = rflg * maxReal;
            return result;
        }
        t = Math.sqrt(rk/z-rk);
        result = rflg*t;
        return result;
    }

    //endregion

    // region Beta, Gamma

    //DotNetHelpers_doc_comment_start  B function for Beta distribution and others 
    //DotNetHelpers_doc_comment_body @param a a parameter 
    //DotNetHelpers_doc_comment_body @param b b parameter 
    //DotNetHelpers_doc_comment_body @return  function value 
    //DotNetHelpers_doc_comment_end
    public final double B(double a, double b)
    {
        return (G(a) * G(b)) / G(a + b);
    }

    //DotNetHelpers_doc_comment_start  vectorial B function for dirichlet distribution 
    //DotNetHelpers_doc_comment_body @param A A vectorial parameter 
    //DotNetHelpers_doc_comment_body @return  the function value 
    //DotNetHelpers_doc_comment_end
    public final double B(List<Double> A)
    {
        double sumA = 0;
        double prodGA = 1;
        for (double a : A)
        {
            sumA += a;
            prodGA *= G(a);
        }
        return prodGA / G(sumA);
    }

    //DotNetHelpers_doc_comment_start  Gamma function for Gamma distribution, Dirichlet distribution and others 
    //DotNetHelpers_doc_comment_body @param k k value (factorial k-1 for integers) 
    //DotNetHelpers_doc_comment_body @return  function value 
    //DotNetHelpers_doc_comment_end
    public final double G(double k)
    {
        return Math.exp(LogG(k));
    }

    private double LogG(double k)
    {
        double x = k, y = k, tmp, ser;
        double[] cof = {76.18009172947146, -86.50532032941677, 24.01409824083091, -1.2317395724500155, 0.1208650973866179e-2, -0.5395239384953e-5};
        tmp = x + 5.5;
        tmp -= (x + 0.5) * Math.log(tmp);
        ser = 1.000000000190015;
        for (int j = 0; j <= 5; j++)
        {
            ser += cof[j] / ++y;
        }
        return -tmp + Math.log(2.5066282746310005 * ser / x);
    }

      //DotNetHelpers_doc_comment_start Returns the incomplete Beta (a, b, x) = Ix(a, b) function. Tested using tables from internet.
    //DotNetHelpers_doc_comment_body @param a a Parameter of the incomplete Beta (a, b, x) function.
    //DotNetHelpers_doc_comment_body @param b b Parameter of the incomplete Beta (a, b, x) function.
    //DotNetHelpers_doc_comment_body @param x x Parameter of the incomplete Beta (a, b, x) function.
    //DotNetHelpers_doc_comment_body @return  calculated value 
      //DotNetHelpers_doc_comment_end
    public final double BetaInc(double a, double b, double x)
    {
        return Stat.betai(a, b, x); //refactor to avoid static call
    }
    //endregion
    
  //region F Distribution

  		public final double FDist(double x, int a, int b)
  		{
  			double result = 0;
  			double w = 0;
  			w = a * x;
  			w = w / (b + w);
  			result = IncompBeta(0.5 * a, 0.5 * b, w);
  			return result;
  		}

  		private double epsilonMachine = 1.11022302462516E-16;

  		public final double IncompBeta(double a, double b, double x)
  		{
  			double result = 0;
  			double t = 0;
  			double xc = 0;
  			double w = 0;
  			double y = 0;
  			int flag = 0;
  			double sg = 0;
  			double big = 0;
  			double biginv = 0;
  			double maxgam = 0;
  			double minlog = 0;
  			double maxlog = 0;

  			big = 4.503599627370496e15;
  			biginv = 2.22044604925031308085e-16;
  			maxgam = 171.624376956302725;
  			minlog = Math.log(epsilonMachine);
  			maxlog = Math.log(Double.MAX_VALUE);
  			if ((double)(x) == (double)(0))
  			{
  				result = 0;
  				return result;
  			}
  			if ((double)(x) == (double)(1))
  			{
  				result = 1;
  				return result;
  			}
  			flag = 0;
  			if ((double)(b * x) <= (double)(1.0) && (double)(x) <= (double)(0.95))
  			{
  				result = IncompBetaps(a, b, x, maxgam);
  				return result;
  			}
  			w = 1.0 - x;
  			if ((double)(x) > (double)(a / (a + b)))
  			{
  				flag = 1;
  				t = a;
  				a = b;
  				b = t;
  				xc = x;
  				x = w;
  			}
  			else
  			{
  				xc = w;
  			}
  			if ((flag == 1 && (double)(b * x) <= (double)(1.0)) && (double)(x) <= (double)(0.95))
  			{
  				t = IncompBetaps(a, b, x, maxgam);
  				if ((double)(t) <= (double)(epsilonMachine))
  				{
  					result = 1.0 - epsilonMachine;
  				}
  				else
  				{
  					result = 1.0 - t;
  				}
  				return result;
  			}
  			y = x * (a + b - 2.0) - (a - 1.0);
  			if ((double)(y) < (double)(0.0))
  			{
  				w = IncompBetafe(a, b, x, big, biginv);
  			}
  			else
  			{
  				w = IncompBetafe2(a, b, x, big, biginv) / xc;
  			}
  			y = a * Math.log(x);
  			t = b * Math.log(xc);
  			if (((double)(a + b) < (double)(maxgam) && (double)(Math.abs(y)) < (double)(maxlog)) && (double)(Math.abs(t)) < (double)(maxlog))
  			{
  				t = Math.pow(xc, b);
  				t = t * Math.pow(x, a);
  				t = t / a;
  				t = t * w;
  				t = t * (Gammafunction(a + b) / (Gammafunction(a) * this.Gammafunction(b)));
  				if (flag == 1)
  				{
  					if ((double)(t) <= (double)(epsilonMachine))
  					{
  						result = 1.0 - epsilonMachine;
  					}
  					else
  					{
  						result = 1.0 - t;
  					}
  				}
  				else
  				{
  					result = t;
  				}
  				return result;
  			}
  			DotNetHelpers.RefObject<Double> tempRef_sg = new DotNetHelpers.RefObject<Double>(sg);
  			DotNetHelpers.RefObject<Double> tempRef_sg2 = new DotNetHelpers.RefObject<Double>(sg);
  			DotNetHelpers.RefObject<Double> tempRef_sg3 = new DotNetHelpers.RefObject<Double>(sg);
  			y = y + t + LnGamma(a + b, tempRef_sg) - this.LnGamma(a, tempRef_sg2) - this.LnGamma(b, tempRef_sg3);
  			sg = tempRef_sg3.argValue;
  			sg = tempRef_sg2.argValue;
  			sg = tempRef_sg.argValue;
  			y = y + Math.log(w / a);
  			if ((double)(y) < (double)(minlog))
  			{
  				t = 0.0;
  			}
  			else
  			{
  				t = Math.exp(y);
  			}
  			if (flag == 1)
  			{
  				if ((double)(t) <= epsilonMachine)
  				{
  					t = 1.0 - epsilonMachine;
  				}
  				else
  				{
  					t = 1.0 - t;
  				}
  			}
  			result = t;
  			return result;
  		}


  		private double IncompBetaps(double a, double b, double x, double maxgam)
  		{
  			double result = 0;
  			double s = 0;
  			double t = 0;
  			double u = 0;
  			double v = 0;
  			double n = 0;
  			double t1 = 0;
  			double z = 0;
  			double ai = 0;
  			double sg = 0;

  			ai = 1.0 / a;
  			u = (1.0 - b) * x;
  			v = u / (a + 1.0);
  			t1 = v;
  			t = u;
  			n = 2.0;
  			s = 0.0;
  			z = epsilonMachine * ai;
  			while ((double)(Math.abs(v)) > (double)(z))
  			{
  				u = (n - b) * x / n;
  				t = t * u;
  				v = t / (a + n);
  				s = s + v;
  				n = n + 1.0;
  			}
  			s = s + t1;
  			s = s + ai;
  			u = a * Math.log(x);
  			if ((double)(a + b) < (double)(maxgam) && (double)(Math.abs(u)) < (double)(Math.log(Double.MAX_VALUE)))
  			{
  				t = Gammafunction(a + b) / (Gammafunction(a) * Gammafunction(b));
  				s = s * t * Math.pow(x, a);
  			}
  			else
  			{
  				DotNetHelpers.RefObject<Double> tempRef_sg = new DotNetHelpers.RefObject<Double>(sg);
  				DotNetHelpers.RefObject<Double> tempRef_sg2 = new DotNetHelpers.RefObject<Double>(sg);
  				DotNetHelpers.RefObject<Double> tempRef_sg3 = new DotNetHelpers.RefObject<Double>(sg);
  				t = LnGamma(a + b, tempRef_sg) - LnGamma(a, tempRef_sg2) - LnGamma(b, tempRef_sg3) + u + Math.log(s);
  				sg = tempRef_sg3.argValue;
  				sg = tempRef_sg2.argValue;
  				sg = tempRef_sg.argValue;
  				if ((double)(t) < (double)(Math.log(epsilonMachine)))
  				{
  					s = 0.0;
  				}
  				else
  				{
  					s = Math.exp(t);
  				}
  			}
  			result = s;
  			return result;
  		}

  		private double IncompBetafe(double a, double b, double x, double big, double biginv)
  		{
  			double result = 0;
  			double xk = 0;
  			double pk = 0;
  			double pkm1 = 0;
  			double pkm2 = 0;
  			double qk = 0;
  			double qkm1 = 0;
  			double qkm2 = 0;
  			double k1 = 0;
  			double k2 = 0;
  			double k3 = 0;
  			double k4 = 0;
  			double k5 = 0;
  			double k6 = 0;
  			double k7 = 0;
  			double k8 = 0;
  			double r = 0;
  			double t = 0;
  			double ans = 0;
  			double thresh = 0;
  			int n = 0;

  			k1 = a;
  			k2 = a + b;
  			k3 = a;
  			k4 = a + 1.0;
  			k5 = 1.0;
  			k6 = b - 1.0;
  			k7 = k4;
  			k8 = a + 2.0;
  			pkm2 = 0.0;
  			qkm2 = 1.0;
  			pkm1 = 1.0;
  			qkm1 = 1.0;
  			ans = 1.0;
  			r = 1.0;
  			n = 0;
  			thresh = 3.0 * epsilonMachine;
  			do
  			{
  				xk = -(x * k1 * k2 / (k3 * k4));
  				pk = pkm1 + pkm2 * xk;
  				qk = qkm1 + qkm2 * xk;
  				pkm2 = pkm1;
  				pkm1 = pk;
  				qkm2 = qkm1;
  				qkm1 = qk;
  				xk = x * k5 * k6 / (k7 * k8);
  				pk = pkm1 + pkm2 * xk;
  				qk = qkm1 + qkm2 * xk;
  				pkm2 = pkm1;
  				pkm1 = pk;
  				qkm2 = qkm1;
  				qkm1 = qk;
  				if ((double)(qk) != (double)(0))
  				{
  					r = pk / qk;
  				}
  				if ((double)(r) != (double)(0))
  				{
  					t = Math.abs((ans - r) / r);
  					ans = r;
  				}
  				else
  				{
  					t = 1.0;
  				}
  				if ((double)(t) < (double)(thresh))
  				{
  					break;
  				}
  				k1 = k1 + 1.0;
  				k2 = k2 + 1.0;
  				k3 = k3 + 2.0;
  				k4 = k4 + 2.0;
  				k5 = k5 + 1.0;
  				k6 = k6 - 1.0;
  				k7 = k7 + 2.0;
  				k8 = k8 + 2.0;
  				if ((double)(Math.abs(qk) + Math.abs(pk)) > (double)(big))
  				{
  					pkm2 = pkm2 * biginv;
  					pkm1 = pkm1 * biginv;
  					qkm2 = qkm2 * biginv;
  					qkm1 = qkm1 * biginv;
  				}
  				if ((double)(Math.abs(qk)) < (double)(biginv) || (double)(Math.abs(pk)) < (double)(biginv))
  				{
  					pkm2 = pkm2 * big;
  					pkm1 = pkm1 * big;
  					qkm2 = qkm2 * big;
  					qkm1 = qkm1 * big;
  				}
  				n = n + 1;
  			} while (n != 300);
  			result = ans;
  			return result;
  		}


  		private double IncompBetafe2(double a, double b, double x, double big, double biginv)
  		{
  			double result = 0;
  			double xk = 0;
  			double pk = 0;
  			double pkm1 = 0;
  			double pkm2 = 0;
  			double qk = 0;
  			double qkm1 = 0;
  			double qkm2 = 0;
  			double k1 = 0;
  			double k2 = 0;
  			double k3 = 0;
  			double k4 = 0;
  			double k5 = 0;
  			double k6 = 0;
  			double k7 = 0;
  			double k8 = 0;
  			double r = 0;
  			double t = 0;
  			double ans = 0;
  			double z = 0;
  			double thresh = 0;
  			int n = 0;

  			k1 = a;
  			k2 = b - 1.0;
  			k3 = a;
  			k4 = a + 1.0;
  			k5 = 1.0;
  			k6 = a + b;
  			k7 = a + 1.0;
  			k8 = a + 2.0;
  			pkm2 = 0.0;
  			qkm2 = 1.0;
  			pkm1 = 1.0;
  			qkm1 = 1.0;
  			z = x / (1.0 - x);
  			ans = 1.0;
  			r = 1.0;
  			n = 0;
  			thresh = 3.0 * epsilonMachine;
  			do
  			{
  				xk = -(z * k1 * k2 / (k3 * k4));
  				pk = pkm1 + pkm2 * xk;
  				qk = qkm1 + qkm2 * xk;
  				pkm2 = pkm1;
  				pkm1 = pk;
  				qkm2 = qkm1;
  				qkm1 = qk;
  				xk = z * k5 * k6 / (k7 * k8);
  				pk = pkm1 + pkm2 * xk;
  				qk = qkm1 + qkm2 * xk;
  				pkm2 = pkm1;
  				pkm1 = pk;
  				qkm2 = qkm1;
  				qkm1 = qk;
  				if ((double)(qk) != (double)(0))
  				{
  					r = pk / qk;
  				}
  				if ((double)(r) != (double)(0))
  				{
  					t = Math.abs((ans - r) / r);
  					ans = r;
  				}
  				else
  				{
  					t = 1.0;
  				}
  				if ((double)(t) < (double)(thresh))
  				{
  					break;
  				}
  				k1 = k1 + 1.0;
  				k2 = k2 - 1.0;
  				k3 = k3 + 2.0;
  				k4 = k4 + 2.0;
  				k5 = k5 + 1.0;
  				k6 = k6 + 1.0;
  				k7 = k7 + 2.0;
  				k8 = k8 + 2.0;
  				if ((double)(Math.abs(qk) + Math.abs(pk)) > (double)(big))
  				{
  					pkm2 = pkm2 * biginv;
  					pkm1 = pkm1 * biginv;
  					qkm2 = qkm2 * biginv;
  					qkm1 = qkm1 * biginv;
  				}
  				if ((double)(Math.abs(qk)) < (double)(biginv) || (double)(Math.abs(pk)) < (double)(biginv))
  				{
  					pkm2 = pkm2 * big;
  					pkm1 = pkm1 * big;
  					qkm2 = qkm2 * big;
  					qkm1 = qkm1 * big;
  				}
  				n = n + 1;
  			} while (n != 300);
  			result = ans;
  			return result;
  		}

  		public final double Gammafunction(double x)
  		{
  			double result = 0;
  			double p = 0;
  			double pp = 0;
  			double q = 0;
  			double qq = 0;
  			double z = 0;
  			int i = 0;
  			double sgngam = 0;

  			sgngam = 1;
  			q = Math.abs(x);
  			if ((double)(q) > (double)(33.0))
  			{
  				if ((double)(x) < (double)(0.0))
  				{
  					p = (int)Math.floor(q);
  					i = (int)Math.round(p);
  					if (i % 2 == 0)
  					{
  						sgngam = -1;
  					}
  					z = q - p;
  					if ((double)(z) > (double)(0.5))
  					{
  						p = p + 1;
  						z = q - p;
  					}
  					z = q * Math.sin(Math.PI * z);
  					z = Math.abs(z);
  					z = Math.PI / (z * Gammastirf(q));
  				}
  				else
  				{
  					z = Gammastirf(x);
  				}
  				result = sgngam * z;
  				return result;
  			}
  			z = 1;
  			while ((double)(x) >= (double)(3))
  			{
  				x = x - 1;
  				z = z * x;
  			}
  			while ((double)(x) < (double)(0))
  			{
  				if ((double)(x) > (double)(-0.000000001))
  				{
  					result = z / ((1 + 0.5772156649015329 * x) * x);
  					return result;
  				}
  				z = z / x;
  				x = x + 1;
  			}
  			while ((double)(x) < (double)(2))
  			{
  				if ((double)(x) < (double)(0.000000001))
  				{
  					result = z / ((1 + 0.5772156649015329 * x) * x);
  					return result;
  				}
  				z = z / x;
  				x = x + 1.0;
  			}
  			if ((double)(x) == (double)(2))
  			{
  				result = z;
  				return result;
  			}
  			x = x - 2.0;
  			pp = 1.60119522476751861407E-4;
  			pp = 1.19135147006586384913E-3 + x * pp;
  			pp = 1.04213797561761569935E-2 + x * pp;
  			pp = 4.76367800457137231464E-2 + x * pp;
  			pp = 2.07448227648435975150E-1 + x * pp;
  			pp = 4.94214826801497100753E-1 + x * pp;
  			pp = 9.99999999999999996796E-1 + x * pp;
  			qq = -2.31581873324120129819E-5;
  			qq = 5.39605580493303397842E-4 + x * qq;
  			qq = -4.45641913851797240494E-3 + x * qq;
  			qq = 1.18139785222060435552E-2 + x * qq;
  			qq = 3.58236398605498653373E-2 + x * qq;
  			qq = -2.34591795718243348568E-1 + x * qq;
  			qq = 7.14304917030273074085E-2 + x * qq;
  			qq = 1.00000000000000000320 + x * qq;
  			result = z * pp / qq;
  			return result;
  		}

  		public final double LnGamma(double x, DotNetHelpers.RefObject<Double> sgngam)
  		{
  			double result = 0;
  			double a = 0;
  			double b = 0;
  			double c = 0;
  			double p = 0;
  			double q = 0;
  			double u = 0;
  			double w = 0;
  			double z = 0;
  			int i = 0;
  			double logpi = 0;
  			double ls2pi = 0;
  			double tmp = 0;

  			sgngam.argValue = 0.0;

  			sgngam.argValue = 1.0;
  			logpi = 1.14472988584940017414;
  			ls2pi = 0.91893853320467274178;
  			if ((double)(x) < (double)(-34.0))
  			{
  				q = -x;
  				DotNetHelpers.RefObject<Double> tempRef_tmp = new DotNetHelpers.RefObject<Double>(tmp);
  				w = LnGamma(q, tempRef_tmp);
  				tmp = tempRef_tmp.argValue;
  				p = (int)Math.floor(q);
  				i = (int)Math.round(p);
  				if (i % 2 == 0)
  				{
  					sgngam.argValue = -1.0;
  				}
  				else
  				{
  					sgngam.argValue = 1.0;
  				}
  				z = q - p;
  				if ((double)(z) > (double)(0.5))
  				{
  					p = p + 1;
  					z = p - q;
  				}
  				z = q * Math.sin(Math.PI * z);
  				result = logpi - Math.log(z) - w;
  				return result;
  			}
  			if ((double)(x) < (double)(13))
  			{
  				z = 1;
  				p = 0;
  				u = x;
  				while ((double)(u) >= (double)(3))
  				{
  					p = p - 1;
  					u = x + p;
  					z = z * u;
  				}
  				while ((double)(u) < (double)(2))
  				{
  					z = z / u;
  					p = p + 1;
  					u = x + p;
  				}
  				if ((double)(z) < (double)(0))
  				{
  					sgngam.argValue = -1.0;
  					z = -z;
  				}
  				else
  				{
  					sgngam.argValue = 1.0;
  				}
  				if ((double)(u) == (double)(2))
  				{
  					result = Math.log(z);
  					return result;
  				}
  				p = p - 2;
  				x = x + p;
  				b = -1378.25152569120859100;
  				b = -38801.6315134637840924 + x * b;
  				b = -331612.992738871184744 + x * b;
  				b = -1162370.97492762307383 + x * b;
  				b = -1721737.00820839662146 + x * b;
  				b = -853555.664245765465627 + x * b;
  				c = 1;
  				c = -351.815701436523470549 + x * c;
  				c = -17064.2106651881159223 + x * c;
  				c = -220528.590553854454839 + x * c;
  				c = -1139334.44367982507207 + x * c;
  				c = -2532523.07177582951285 + x * c;
  				c = -2018891.41433532773231 + x * c;
  				p = x * b / c;
  				result = Math.log(z) + p;
  				return result;
  			}
  			q = (x - 0.5) * Math.log(x) - x + ls2pi;
  			if ((double)(x) > (double)(100000000))
  			{
  				result = q;
  				return result;
  			}
  			p = 1 / (x * x);
  			if ((double)(x) >= (double)(1000.0))
  			{
  				q = q + ((7.9365079365079365079365 * 0.0001 * p - 2.7777777777777777777778 * 0.001) * p + 0.0833333333333333333333) / x;
  			}
  			else
  			{
  				a = 8.11614167470508450300 * 0.0001;
  				a = -(5.95061904284301438324 * 0.0001) + p * a;
  				a = 7.93650340457716943945 * 0.0001 + p * a;
  				a = -(2.77777777730099687205 * 0.001) + p * a;
  				a = 8.33333333333331927722 * 0.01 + p * a;
  				q = q + a / x;
  			}
  			result = q;
  			return result;
  		}

  		private double Gammastirf(double x)
  		{
  			double result = 0;
  			double y = 0;
  			double w = 0;
  			double v = 0;
  			double stir = 0;

  			w = 1 / x;
  			stir = 7.87311395793093628397E-4;
  			stir = -2.29549961613378126380E-4 + w * stir;
  			stir = -2.68132617805781232825E-3 + w * stir;
  			stir = 3.47222221605458667310E-3 + w * stir;
  			stir = 8.33333333333482257126E-2 + w * stir;
  			w = 1 + w * stir;
  			y = Math.exp(x);
  			if ((double)(x) > (double)(143.01608))
  			{
  				v = Math.pow(x, 0.5 * x - 0.25);
  				y = v * (v / y);
  			}
  			else
  			{
  				y = Math.pow(x, x - 0.5) / y;
  			}
  			result = 2.50662827463100050242 * y * w;
  			return result;
  		}

  		//endregion

  		//endregion

  		//region Homoscedasticity

  		/**  Log returns calculation 
  		 @param serie time series 
  		 @param abs if absolute log returns or not 
  		 @return  log returns 
  		*/
  		public final ArrayList<Double> LogReturns(List<Double> serie, boolean abs)
  		{
  			double min = Double.MAX_VALUE;
  			for (double val : serie)
  			{
  				if (val < min)
  				{
  					min = val;
  				}
  			}
  			if (min <= 0)
  			{
  				for (int i = 1;i < serie.size();i++)
  				{
  					serie.set(i, serie.get(i) - min + 1);
  				}
  			}
  			ArrayList<Double> logReturns = new ArrayList<Double>();
  			for (int i = 1;i < serie.size();i++)
  			{
  				if (serie.get(i) < 0)
  				{
  					serie.set(i, 0.0);
  				}
  				if (abs)
  				{
  					logReturns.add(Math.abs(Math.log(serie.get(i)) - Math.log(serie.get(i - 1))));
  				}
  				else
  				{
  					logReturns.add(Math.log(serie.get(i)) - Math.log(serie.get(i - 1)));
  				}
  			}
  			return logReturns;
  		}

  		public final double WeightedMean(List<Double> values, List<Double> weigths)
  		{
  			if (values.size() != weigths.size())
  			{
  				throw new RuntimeException("Error. Values and weights must have same count");
  			}
  			double sum = 0;
  			for (int i = 0; i < values.size(); i++)
  			{
  				sum += values.get(i) * weigths.get(i);
  			}
  			return sum / (double)values.size();
  		}

  		public final double WeightedMean(List<Double> values, int firstIndex, int lastIndex, List<Double> weights)
  		{
  			if (lastIndex - firstIndex + 1 != weights.size())
  			{
  				throw new RuntimeException("Error. Values and weights must have same count");
  			}
  			if (firstIndex < 0 || lastIndex >= values.size())
  			{
  				throw new RuntimeException("Error. intex out of range");
  			}
  			double sum = 0;
  			double wCount = 0;
  			for (int i = firstIndex; i <= lastIndex; i++)
  			{
  				sum += values.get(i) * weights.get(i - firstIndex);
  				wCount += weights.get(i - firstIndex);
  			}
  			return sum / wCount;
  		}

  		public final ArrayList<Double> GetInversePerc(List<Double> percs)
  		{
  			ArrayList<Double> inversePerc = new ArrayList<Double>();
  			for (int i = 0; i < percs.size(); i++)
  			{
  				inversePerc.add(1 - percs.get(i));
  			}
  			return inversePerc;
  		}
  		/**  Weighted variance 
  		 @param serie serie for calculation 
  		 @param s standard deviation 
  		 @return  weighted variance list 
  		*/
  		public final ArrayList<Double> WeightedVar(List<Double> serie, double s)
  		{
  			ArrayList<Double> wv = new ArrayList<Double>();
  			for (int i = 0;i < serie.size();i++)
  			{
  				wv.add(WeightedVar(serie, i, s));
  			}
  			return wv;
  		}

  		/**  Weighted variance from a certain index 
  		 @param serie serie for calculation 
  		 @param index first index 
  		 @param s standard deviation 
  		 @return  weighted variance list 
  		*/
  		public final double WeightedVar(List<Double> serie, int index, double s)
  		{
  			NormalDist nd = new NormalDist();
  			ArrayList<Double> dist = new ArrayList<Double>();
  			double eps = 0.01;
  			double w = 1;
  			int m, x = 0;
  			m = index;
  			x = 0;
  			while (w > eps && m - x >= 0)
  			{
  				w = nd.pNormDeriv(m - x, m, s);
  				for (int i = 0;i < w * 100;i++)
  				{
  					if (serie.get(m - x) < 0)
  					{
  						dist.add(0.0);
  					}
  					else
  					{
  						dist.add(serie.get(m - x));
  					}
  				}
  				x++;
  			}
  			x = 0;
  			while (w > eps && m + x < serie.size())
  			{
  				w = nd.pNormDeriv(m + x, m, s);
  				for (int i = 0;i < w * 100;i++)
  				{
  					if (serie.get(m + x) < 0)
  					{
  						dist.add(0.0);
  					}
  					else
  					{
  						dist.add(serie.get(m + x));
  					}
  				}
  				x++;
  			}
  			return Variance(dist);
  		}

  		//endregion

  		//region ACF

  		/**  Autocovariance calculation 
  		 @param serie time series 
  		 @param lag lag for autocovariance 
  		 @return  autocovariance value C_h 
  		*/
  		public final double Autocov(List<Double> serie, int lag)
  		{
  			double m = Mean(serie);
  			double n = serie.size();
  			double ac = 0;
  			for (int i = 0;i < serie.size() - lag;i++)
  			{
  				ac += (serie.get(i) - m) * (serie.get(i + lag) - m);
  			}
  			return ac / n;
  		}


  		/**  ACF calculation 
  		 @param serie time series 
  		 @param lag lag for autocorrelation 
  		 @return  autocorrelation coefficient R_h 
  		*/
  		public final double Autocorr(List<Double> serie, int lag)
  		{
  			double m = Mean(serie);
  			double ac = 0;
  			double va = 0;
  			double n = serie.size();
  			for (int i = 0;i < n;i++)
  			{
  				if (i < n - lag)
  				{
  					ac += (serie.get(i) - m) * (serie.get(i + lag) - m);
  				}
  				va += Math.pow((serie.get(i) - m),2);
  			}
  			return ac / va;
  		}

  		/**  Sample autocorrelation function series 
  		 @param serie time series 
  		 @param maxLag lag for autocorrelation (1 to lag) 
  		 @return  autocorrelation time series 
  		*/
  		public final ArrayList<Double> ACF(List<Double> serie, int maxLag)
  		{
  			ArrayList<Double> autocorrelation = new ArrayList<Double>();
  			for (int l = 0;l <= maxLag;l++)
  			{
  				autocorrelation.add(Autocorr(serie, l));
  			}
  			return autocorrelation;
  		}

  		public final int CalculateMaxAutocorrLag(ArrayList<Double> serie, double minAutoCorr, int minLag, int maxLag, double epsilon)
  		{
  			double bestAc = -Double.MAX_VALUE;
  			int bestLag = -1;
  			double ac;
  			for (int lag = minLag; lag <= maxLag; lag++)
  			{
  				ac = Autocorr(serie, lag);
  				if (ac > bestAc + epsilon)
  				{
  					bestAc = ac;
  					bestLag = lag;
  				}
  			}
  			if (bestAc > minAutoCorr)
  			{
  				return bestLag;
  			}
  			else
  			{
  				return -1;
  			}
  		}

  		public final ArrayList<Integer> CalculateEnoughAutocorrLags(ArrayList<Double> serie, double minAutoCorr, int minLag, int maxLag)
  		{
  			ArrayList<Integer> enoughLags = new ArrayList<Integer>();
  			double ac;
  			for (int lag = minLag; lag <= maxLag; lag++)
  			{
  				ac = Autocorr(serie, lag);
  				if (ac >= minAutoCorr)
  				{
  					enoughLags.add(lag);
  				}
  			}
  			return enoughLags;
  		}


  		public final ArrayList<Double> CalculateLagSerie(ArrayList<Double> serie, int lag, int mw)
  		{
  			ArrayList<Double> lagSerie = new ArrayList<Double>();
  			for (int i = serie.size() - lag; i > 0; i = i - lag)
  			{
  				for (int j = i - mw / 2; j < i + (mw - mw / 2); j++)
  				{
  					lagSerie.add(serie.get(j));
  				}
  			}
  			return lagSerie;
  		}

  		//endregion

  		//region Partial autocorrelation

  		/**  Calculates the sample Partial autocorrelation function using Durbin's algorithm 
  		 @param acfs autocorrelation function series 
  		 @return  partial autocorrelation function series 
  		*/
  		public final ArrayList<Double> PACF(List<Double> acfs) 	{
  			double[] pacfs = null;
  			double[] var = null;
  			double[][] alpha = null;

  			// get acfs
  			int lags = acfs.size();
  			pacfs = new double[lags];
  			var = new double[lags];
  			alpha = new double[lags][];
  			for (int i = 0;i < lags;i++)
  			{
  				alpha[i] = new double[lags];
  			}

  			//values 0-1
  			pacfs[0] = acfs.get(0);
  			var[0] = 0;
  			alpha[0][0] = 0;
  			alpha[1][0] = 0;
  			alpha[0][1] = 0;
  			pacfs[1] = acfs.get(1);
  			alpha[1][1] = pacfs[1];
  			var[1] = 1 - Math.pow(pacfs[1], 2);

  			//values 2-lags with iterative formula
  			for (int k = 1;k < lags - 1;k++)
  			{
  				double sum1 = 0;
  				for (int j = 1;j <= k;j++)
  				{
  					sum1 += alpha[j][k] * acfs.get(k + 1 - j);
  				}

  				pacfs[k + 1] = (acfs.get(k + 1) - sum1) / var[k];
  				alpha[k + 1][k + 1] = pacfs[k + 1];

  				//now determine remaining alphas
  				for (int j = 1;j <= k;j++)
  				{
  					alpha[j][k + 1] = alpha[j][k] - pacfs[k + 1] * alpha[k + 1 - j][k];
  				}
  				var[k + 1] = var[k] * (1 - Math.pow(pacfs[k + 1], 2));
  			}
  			ArrayList<Double> arrL = new ArrayList<Double>();
  			for(int i=0;i<arrL.size();i++) { arrL.add(pacfs[i]); }
  			return(arrL);
  		}




  		//endregion

  		//region Standard Score

  		public final ArrayList<Double> StandardScore(List<Double> values)
  		{
  			ArrayList<Double> ss = new ArrayList<Double>();
  			double mean = Mean(values);
  			double sd = StDev(values);
  			for (int i = 0; i < values.size(); i++)
  			{
  				ss.add((values.get(i) - mean) / sd);
  			}
  			return ss;
  		}

  		public final ArrayList<Double> StandardScore(List<Double> values, int firstIndex, int lastIndex)
  		{
  			if (firstIndex < 0 || lastIndex >= values.size() || lastIndex < firstIndex)
  			{
  				throw new RuntimeException("Error. intex out of range");
  			}
  			ArrayList<Double> ss = new ArrayList<Double>();
  			if (firstIndex == lastIndex)
  			{
  				ss.add(1.0);
  				return ss;
  			}
  			double mean = Mean(values);
  			double sd = StDev(values);
  			for (int i = firstIndex; i <= lastIndex; i++)
  			{
  				ss.add((values.get(i) - mean) / sd);
  			}
  			return ss;
  		}

  		//endregion
  		//endregion

}
