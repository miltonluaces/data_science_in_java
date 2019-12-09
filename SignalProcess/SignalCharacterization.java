package SignalProcess;

import java.util.*;

//region Imports

import SignalProcess.*;
import NumCalc.*;
import ClassicalStat.*;
import Clustering.*;

//endregion


//Tangible_doc_comment_start Class for signal characterization: trend removal, noise filtering, outlier filtering, stationarity, statistics 
//Tangible_doc_comment_end
public class SignalCharacterization
{

    //region Fields

    //calculation classes
    private Wavelets wav;
    private Polynomic poly;
    private Statistics stat;
    private Functions func;
    private NormalDist ns;
    private Derivatives drv;
    private HomogeneityClustering hc;
    private ArrayList<HomogeneityClustering.StatCluster> clusters;
    private HomogeneityClustering.StatCluster firstCluster;
    private HypoTest hypoTest;
    private double[] weights;

    //parameters
    private double varCoeffThreshold;
    private boolean filterOutliers;
    private int minFreqForOutlierFiltering;
    private int minCount;
    private int minCountForClustering;
    private double outlierThreshold;
    private double dailyOutlierThreshold;
    private double maxNorm;
    private double minReg;
    private double maxReg;
    private int forgetInitPeriods;
    private int forgetEndPeriods;
    private double forgetEndProportion;
    private double minLeadTimeObsProportion;
    private double pValue = 0.05;
    private double fc = 0;
    private int minSpan;
    private double loadFactorThreshold;

    //collections
    private ArrayList<Double> origSerie;
    private ArrayList<Double> filteredSerie;
    private ArrayList<Double> serie;
    //private List<double> serieDrv1;
    private ArrayList<Double> serieDrv2;
    private ArrayList<Double> regres;
    private ArrayList<Double> bias;

    private ArrayList<Double> normSerie;
    private ArrayList<Double> normRegres;

    private ArrayList<Double> means;
    private ArrayList<Double> stDevs;

    private TreeMap<Integer, Double> freqsDict;

    //results
    private double min;
    private double max;
    private double minOutlierFiltered;
    private double maxValueNotFiltered;
    private ArrayList<Integer> outlierIndexes;
    private int firstIndexToFcst;
    private ArrayList<ResultType> results;
    private double statPeriodMean;
    private double minNonZeroFreq;
    private double significance;
    private HomogeneityClustering.TestType testType;

    private double minFc;
    private boolean sparse;
    private int lag;
    private boolean continuous;
    private int span;

    private double minWall;
    private int nSpans;

    private double rootExpForClusterWeighting;
    private double ratioMvThreshold;

    private boolean trace = false;

    //endregion

    //region Constructor

    //Tangible_doc_comment_start  Constructor 
    //Tangible_doc_comment_body @param minCount minimum quantity of values 
    //Tangible_doc_comment_body @param varCoeffThreshold variance coefficent threshold 
    //Tangible_doc_comment_body @param filterOutliers if outliers should be filtered 
    //Tangible_doc_comment_body @param minFreqForOutlierFiltering minimum frequency to allow outlier filtering 
    //Tangible_doc_comment_body @param outlierThreshold Threshold for oultier filtering 
    //Tangible_doc_comment_body @param maxNorm Maximum normalized value (i.e. 100) 
    //Tangible_doc_comment_body @param maxDiff Maximum mean difference allowed between two clusters 
    //Tangible_doc_comment_body @param maxStDev Maximum standard deviation allowed in a cluster (after join) 
    //Tangible_doc_comment_body @param pValue pValue por hypothesis tests 
    //Tangible_doc_comment_end
    public SignalCharacterization(int minCount, int minCountForClustering, double varCoeffThreshold, boolean filterOutliers, int minFreqForOutlierFiltering, double outlierThreshold, double dailyOutlierThreshold, double maxNorm, double maxDiff, double maxStDev, double pValue)
    {

        this.minCount = minCount;
        this.minCountForClustering = minCountForClustering;
        this.varCoeffThreshold = varCoeffThreshold;
        this.filterOutliers = filterOutliers;
        this.minFreqForOutlierFiltering = minFreqForOutlierFiltering;
        this.outlierThreshold = outlierThreshold;
        this.dailyOutlierThreshold = dailyOutlierThreshold;
        this.maxNorm = maxNorm;

        this.freqsDict = new TreeMap<Integer, Double>();
        this.results = new ArrayList<ResultType>();

        this.poly = new Polynomic();
        this.stat = new Statistics();
        this.func = new Functions();
        this.wav = new Wavelets();
        this.ns = new NormalDist();
        this.drv = new Derivatives();
        this.hypoTest = new HypoTest(0.05, 5);
        this.testType = HomogeneityClustering.TestType.Statistics;
        //this.testType = HomogeneityClustering.TestType.Wilcoxon;
        this.hc = new HomogeneityClustering(maxDiff, maxStDev, minCountForClustering, pValue, testType);
        this.statPeriodMean = -1.0;
        this.clusters = new ArrayList<HomogeneityClustering.StatCluster>();
        this.minNonZeroFreq = Double.MAX_VALUE;
        this.pValue = pValue;
        this.minSpan = 7;
        this.minFc = 0.8;
        this.minCount = 5;
        this.minCountForClustering = 5;
        this.lag = -1;
        this.span = 1;
        this.minWall = 2;
        this.rootExpForClusterWeighting = 1;
        this.ratioMvThreshold = 0.2;
    }

    //endregion

    //region Properties

    //Tangible_doc_comment_start  Minimum quantity of values allowed 
    //Tangible_doc_comment_end
    public final int getMinCount()
    {
        return minCount;
    }
    public final void setMinCount(int value)
    {
        minCount = value;
    }

    public final int getMinCountForClustering()
    {
        return minCountForClustering;
    }
    public final void setMinCountForClustering(int value)
    {
        minCountForClustering = value;
        hc.MinCount = value;
    }

    //Tangible_doc_comment_start  load factor of grouped non-zero values on total grouped values 
    //Tangible_doc_comment_end
    public final double getLoadFactor()
    {
        return wav.LoadFactor;
    }

    //Tangible_doc_comment_start  Minimum value 
    //Tangible_doc_comment_end
    public final double getMin()
    {
        return min;
    }

    //Tangible_doc_comment_start  Maximum value 
    //Tangible_doc_comment_end
    public final double getMax()
    {
        return max;
    }

    //Tangible_doc_comment_start  Range of values 
    //Tangible_doc_comment_end
    public final double getRange()
    {
        return max-min;
    }

    //Tangible_doc_comment_start  Min value of a filtered outlier 
    //Tangible_doc_comment_end
    public final double getMinOutlierFiltered()
    {
        return minOutlierFiltered;
    }

    //Tangible_doc_comment_start  Max value of a non-filtered value 
    //Tangible_doc_comment_end
    public final double getMaxValueNotFiltered()
    {
        return maxValueNotFiltered;
    }

    //Tangible_doc_comment_start  Filtered time series 
    //Tangible_doc_comment_end
    public final ArrayList<Double> getFilteredSerie()
    {
        return filteredSerie;
    }

    ///// <summary> First derivative time series </summary>
    //public List<double> SerieDrv1 { 
        //get { return serieDrv1; } 
    //}

    //Tangible_doc_comment_start  Second derivative time series 
    //Tangible_doc_comment_end
    public final ArrayList<Double> getSerieDrv2()
    {
        return serieDrv2;
    }

    //Tangible_doc_comment_start  Index of first period to forecast 
    //Tangible_doc_comment_end
    public final int getFirstIndexToFcst()
    {
        return firstIndexToFcst;
    }

    //Tangible_doc_comment_start  Original time series 
    //Tangible_doc_comment_end
    public final ArrayList<Double> getOrigSerie()
    {
        return origSerie;
    }
    public final void setOrigSerie(ArrayList<Double> value)
    {
        origSerie = value;
    }

    //Tangible_doc_comment_start  Current time series 
    //Tangible_doc_comment_end
    public final ArrayList<Double> getSerie()
    {
        return serie;
    }
    public final void setSerie(ArrayList<Double> value)
    {
        serie = value;
    }

    //Tangible_doc_comment_start  Regression time series 
    //Tangible_doc_comment_end
    public final ArrayList<Double> getRegres()
    {
        return regres;
    }

    //Tangible_doc_comment_start  Bias time series 
    //Tangible_doc_comment_end
    public final ArrayList<Double> getBias()
    {
        return bias;
    }

    //Tangible_doc_comment_start  Indexes of filtered outliers 
    //Tangible_doc_comment_end
    public final ArrayList<Integer> getOutlierIndexes()
    {
        return outlierIndexes;
    }

    //Tangible_doc_comment_start  List of results (warinings or errors) 
    //Tangible_doc_comment_end
    public final ArrayList<ResultType> getResults()
    {
        return results;
    }

    //Tangible_doc_comment_start  Data frequencies in a dictionary 
    //Tangible_doc_comment_end
    public final TreeMap<Integer, Double> getFreqsDict()
    {
        return freqsDict;
    }

    //Tangible_doc_comment_start  Threshold for outlier filtering 
    //Tangible_doc_comment_end
    public final double getOutlierThreshold()
    {
        return outlierThreshold;
    }
    public final void setOutlierThreshold(double value)
    {
        outlierThreshold = value;
    }

    //Tangible_doc_comment_start  Threshold for daily outlier filtering 
    //Tangible_doc_comment_end
    public final double getDailyOutlierThreshold()
    {
        return dailyOutlierThreshold;
    }
    public final void setDailyOutlierThreshold(double value)
    {
        dailyOutlierThreshold = value;
    }

    //Tangible_doc_comment_start  Maximum mean difference allowed between two clusters 
    //Tangible_doc_comment_end
    public final double getMaxDiff()
    {
        return hc.MaxDiff;
    }
    public final void setMaxDiff(double value)
    {
        hc.MaxDiff = value;
    }

    //Tangible_doc_comment_start  Maximum standard deviation difference allowed between two clusters 
    //Tangible_doc_comment_end
    public final double getMaxStDevDiff()
    {
        return hc.MaxStDevDiff;
    }
    public final void setMaxStDevDiff(double value)
    {
        hc.MaxStDevDiff = value;
    }

    //Tangible_doc_comment_start  Maximum difference allowed within a cluster (after join) 
    //Tangible_doc_comment_end
    public final double getMaxStDev()
    {
        return hc.MaxStDev;
    }
    public final void setMaxStDev(double value)
    {
        hc.MaxStDev = value;
    }

    //Tangible_doc_comment_start  List of cluster means 
    //Tangible_doc_comment_end
    public final ArrayList<Double> getMeans()
    {
        return means;
    }

    //Tangible_doc_comment_start  List of cluster standard deviations 
    //Tangible_doc_comment_end
    public final ArrayList<Double> getStDevs()
    {
        return stDevs;
    }

    //Tangible_doc_comment_start  Array of weights for clusters 
    //Tangible_doc_comment_end
    public final double[] getWeights()
    {
        return weights;
    }
    public final void setWeights(double[] value)
    {
        weights = value;
    }

    //Tangible_doc_comment_start  Mean of stationary period 
    //Tangible_doc_comment_end
    public final double getStatPeriodMean()
    {
        return statPeriodMean;
    }

    //Tangible_doc_comment_start  Minimum leadtime observations in proportion with total observations 
    //Tangible_doc_comment_end
    public final double getMinLeadTimeObsProportion()
    {
        return minLeadTimeObsProportion;
    }
    public final void setMinLeadTimeObsProportion(double value)
    {
        minLeadTimeObsProportion = value;
    }

    //Tangible_doc_comment_start  List of clusters 
    //Tangible_doc_comment_end
    public final ArrayList<HomogeneityClustering.StatCluster> getClusters()
    {
        return clusters;
    }

    //Tangible_doc_comment_start  First index of stationary period 
    //Tangible_doc_comment_end
    public final int getStatPeriodFirstIndex()
    {
        return firstIndexToFcst;
    }
    public final void setStatPeriodFirstIndex(int value)
    {
        firstIndexToFcst = value;
    }

    //Tangible_doc_comment_start  Start of forgetting period (from the present to the past) 
    //Tangible_doc_comment_end
    public final int getForgetInitPeriods()
    {
        return forgetInitPeriods;
    }
    public final void setForgetInitPeriods(int value)
    {
        forgetInitPeriods = value;
    }

    //Tangible_doc_comment_start  End of forgetting period (from the present to the past)  
    //Tangible_doc_comment_end
    public final int getForgetEndPeriods()
    {
        return forgetEndPeriods;
    }
    public final void setForgetEndPeriods(int value)
    {
        forgetEndPeriods = value;
    }

    //Tangible_doc_comment_start  Forgetting proportion beyond end of stationary period (from the present to the past) 
    //Tangible_doc_comment_end
    public final double getForgetEndProportion()
    {
        return forgetEndProportion;
    }
    public final void setForgetEndProportion(double value)
    {
        forgetEndProportion = value;
    }

    //Tangible_doc_comment_start  Number of elements of first cluster 
    //Tangible_doc_comment_end
    public final int getFirstClusterCount()
    {
        return firstCluster.size();
    }

    //Tangible_doc_comment_start  Type of statistical test for joins 
    //Tangible_doc_comment_end
    public final HomogeneityClustering.TestType getTest()
    {
        return testType;
    }
    public final void setTest(HomogeneityClustering.TestType value)
    {
        testType = value;
        hc.Test = value;
    }

    //Tangible_doc_comment_start  Minimum Load Factor 
    //Tangible_doc_comment_end
    public final double getMinFc()
    {
        return minFc;
    }
    public final void setMinFc(double value)
    {
        minFc = value;
    }

    //Tangible_doc_comment_start  Cut-off for trimming 
    //Tangible_doc_comment_end
    public final double getCutOff()
    {
        return hc.CutOff;
    }
    public final void setCutOff(double value)
    {
        hc.CutOff = value;
    }

    //Tangible_doc_comment_start  load factor threshold 
    //Tangible_doc_comment_end
    public final double getLoadFactorThreshold()
    {
        return loadFactorThreshold;
    }
    public final void setLoadFactorThreshold(double value)
    {
        loadFactorThreshold = value;
    }

    //Tangible_doc_comment_start  if it is sparse 
    //Tangible_doc_comment_end
    public final boolean getSparse()
    {
        return sparse;
    }
    public final void setSparse(boolean value)
    {
        sparse = value;
    }

    //Tangible_doc_comment_start  lag for seasonality 
    //Tangible_doc_comment_end
    public final int getLag()
    {
        return lag;
    }
    public final void setLag(int value)
    {
        lag = value;
    }

    //Tangible_doc_comment_start  Exponent for non linear weighting of clusters (1-25) 
    //Tangible_doc_comment_end
    public final double getRootExpForClusterWeighting()
    {
        return rootExpForClusterWeighting;
    }
    public final void setRootExpForClusterWeighting(double value)
    {
        rootExpForClusterWeighting = value;
    }

    //Tangible_doc_comment_start  Ratio mean/variance for sparse time series (threshold) 
    //Tangible_doc_comment_end
    public final double getRatioMvThreshold()
    {
        return ratioMvThreshold;
    }
    public final void setRatioMvThreshold(double value)
    {
        ratioMvThreshold = value;
    }

    //endregion

    //region Public Methods

    //Tangible_doc_comment_start  Main calculation method 
    //Tangible_doc_comment_body @param minCount minimum number of cases 
    //Tangible_doc_comment_end
    public final void Characterization(int minCount)
    {
        SetRegressionAndBias();
        FilterOutliers();
        NormalizeCollections();
        CalculateWeights();
        HomogeneityClusteringTest(firstIndexToFcst);
        SetClusterWeights();
    }

    //Tangible_doc_comment_start  Calculate load factor for sparse calculations 
    //Tangible_doc_comment_body @param origSerie original serie (not grouped) 
    //Tangible_doc_comment_end
    public final void CalculateLoadFactor(ArrayList<Double> origSerie)
    {
        int nonZerosCount = 0;
        for (double val : origSerie)
        {
            if (val > 0)
            {
                nonZerosCount++;
            }
        }
        fc = (double)nonZerosCount / (double)origSerie.size();
    }

    //Tangible_doc_comment_start  Calculate index of first homogeneous period 
    //Tangible_doc_comment_body @return  the index 
    //Tangible_doc_comment_end
    public final int CalculateHomogeneousInitialIndex(double alpha, double highAlpha)
    {
        int maxLag = 12;
        int initialIndex = 0;
        if(clusters.size() < 2)
        {
            return initialIndex;
        }

        double pValue1 = -1;
        double pValue2 = -1;
        for(int i=clusters.size()-1;i>0;i--)
        {
            initialIndex = clusters.get(i).FirstIndex;
            tangible.RefObject<Double> tempRef_pValue1 = new tangible.RefObject<Double>(pValue1);
            tangible.RefObject<Double> tempRef_pValue2 = new tangible.RefObject<Double>(pValue2);
            if(clusters.get(i-1).size() > 1 && clusters.get(i).size() > 1 && !(hc.AreHomogeneous(this.origSerie, clusters.get(i).FirstIndex, tempRef_pValue1) || hc.AreHomogeneous(origSerie, clusters.get(i-1).FirstIndex, clusters.get(i-1).LastIndex, clusters.get(i).FirstIndex, clusters.get(i).LastIndex, tempRef_pValue2)) && origSerie.size() - initialIndex >= minCountForClustering)
            {
            pValue2 = tempRef_pValue2.argValue;
            pValue1 = tempRef_pValue1.argValue;
                break;
            }
        else
        {
            pValue2 = tempRef_pValue2.argValue;
            pValue1 = tempRef_pValue1.argValue;
        }
        }
        if (initialIndex == 0 || pValue1 < highAlpha || pValue2 < highAlpha)
        {
            return initialIndex;
        }

        ArrayList<Double> acf = stat.ACF(serie, maxLag);
        int p = hypoTest.TestMovingAverage(acf, alpha);
        if (p != -1)
        {
            return 0;
        }

        return initialIndex;
    }


    //Tangible_doc_comment_start  Calculate index of first homogeneous period 
    //Tangible_doc_comment_body @return  the index 
    //Tangible_doc_comment_end
    @Deprecated
    public final int CalculateHomogeneousInitialIndexObsol()
    {
        int initialIndex = 0;
        if(clusters.size() < 2)
        {
            return initialIndex;
        }
        HomogeneityClustering.StatCluster statCluster = null;
        HomogeneityClustering.StatCluster restCluster;
        initialIndex = clusters.get(clusters.size()-1).FirstIndex;
        for(int i=clusters.size()-1;i>0;i--)
        {
            if(statCluster == null)
            {
                statCluster = clusters.get(i);
            }
            else
            {
                statCluster.Join(clusters.get(i));
            }

            restCluster = new HomogeneityClustering.StatCluster(clusters.get(i-1));
            //if(i>=2) {  for(int j=i-2;j>=0;j--) { restCluster.Join(clusters[j]); } }

            if(!hc.AreHomogeneous(clusters.get(i), restCluster) || serie.size() - statCluster.FirstIndex >= minCountForClustering)
            {
                break;
            }
            initialIndex = statCluster.FirstIndex;
        }
        return initialIndex;
    }

    @Deprecated
    public final int CalculateHomogeneousInitialIndexObsol2()
    {
        int initialIndex = 0;
        if(clusters.size() < 2)
        {
            return initialIndex;
        }
        HomogeneityClustering.StatCluster restCluster;
        initialIndex = clusters.get(0).FirstIndex;
        for(int i=0;i<clusters.size()-1;i++)
        {
            restCluster = new HomogeneityClustering.StatCluster(clusters.get(i+1));
            for(int j=i+2;j<clusters.size();j++)
            {
                restCluster.Join(clusters.get(j));
            }
            if(hc.AreHomogeneous(clusters.get(i), restCluster) || serie.size() - clusters.get(i).size() < minCountForClustering)
            {
                break;
            }
            initialIndex = restCluster.FirstIndex;
        }
        return initialIndex;
    }

    public final void JoinSmallClusters(int minSize, double pValueThreshold)
    {
        for(int i=0;i<clusters.size();i++)
        {
            if(clusters.get(i).size() > minSize)
            {
                continue;
            }
            double pVal1 = 0;
            double pVal2 = 0;
            if(i>0)
            {
                pVal1 = hc.TestHomogeneity(clusters.get(i-1), clusters.get(i));
            }
            if(i<clusters.size()-1)
            {
                pVal2 = hc.TestHomogeneity(clusters.get(i), clusters.get(i+1));
            }
            if(pVal1 >= pValueThreshold && pVal1 > pVal2)
            {
                clusters.get(i).Join(clusters.get(i-1));
            }
            else if(pVal2 >= pValueThreshold && pVal2 > pVal1)
            {
                clusters.get(i).Join(clusters.get(i+1));
            }
        }
    }

    public final double TestHomogeneity(HomogeneityClustering.StatCluster sc1, HomogeneityClustering.StatCluster sc2)
    {
        return hc.TestHomogeneity(sc1, sc2);
    }

    //endregion

    //region Private Methods

    //region Load Data

    public final void LoadData(ArrayList<Double> origSerie, int firstIndex, int span, int nSpans, int minCountForClustering, boolean continuous)
    {
        this.continuous = continuous;
        this.nSpans = nSpans;
        this.span =+ span;
        if (minCountForClustering != -1)
        {
            this.minCountForClustering = minCountForClustering;
        }
        if (minCountForClustering < 3)
        {
            minCountForClustering = (int)Math.min(origSerie.size(), 3);
        }
        if (continuous)
        {
            forgetInitPeriods = forgetInitPeriods * span;
            forgetEndPeriods = forgetEndPeriods * span;
            this.minCountForClustering = minCountForClustering * span;
            this.minWall = minWall * span;
        }
        this.results.clear();
        this.origSerie = origSerie;
        FilterDailyOutliers(origSerie, firstIndex);

        clusters.clear();
        wav.Group(this.origSerie, 0, span, continuous);
        this.serie = wav.Data;

        //is sparse? 
        //double lf = GetLoadFactor(serie);
        double lf = GetWeightedLoadFactor(serie);
        if (this.getLoadFactorThreshold() == 0 || (lf > 0.05 && (lf >= this.getLoadFactorThreshold() || clusters.size() > 1)))
        {
            sparse = false;
        }
        else
        {
            double m = stat.Mean(serie);
            double v = stat.Variance(serie);
            if (v > 0 && m / v > ratioMvThreshold)
            {
                sparse = false;
            }
            else
            {
                sparse = true;
            }
        }
        SetMinMaxAndFreqs(serie);

        if (sparse)
        {
            return;
        }

        int maxSpan = (int)((double)origSerie.size() / 4.00); //max span to warrant 4 grouped data
        if (span > maxSpan && maxSpan >= 7)
        {
            span = maxSpan;
        }
    }

    private ArrayList<Double> Group(ArrayList<Double> origSerie, int span, int firstIndex)
    {
        ArrayList<Double> serie = new ArrayList<Double>();
        double groupedValue = 0;

        for (int i = firstIndex; i < origSerie.size(); i++)
        {
            groupedValue += origSerie.get(i);
            if ((i + 1) % span == 0)
            {
                serie.add(groupedValue);
                groupedValue = 0;
            }
        }
        if (origSerie.size() % span != 0)
        {
            serie.add(groupedValue);
        }
        return serie;
    }

    //endregion

    //region Regression and Bias

    //Tangible_doc_comment_start  Set wavelet regression and calculate bias 
    //Tangible_doc_comment_end
    public final void SetRegressionAndBias()
    {
        if (freqsDict.size() == 1 || serie.size() <= 2)
        {
            return;
        }
        else if(serie.size() <= 4)
        {
            regres = new ArrayList<Double>();
            regres.addAll(serie);
        }
        else
        {
            wav.CalcFreqDescomposition();
            regres = wav.GetFreqDescomposition(wav.Coeffs.size()-2);
        }
        bias = func.Substract(serie, regres, false, false);
    }

    //endregion

    //region Outlier Filtering

    //region Grouped filtering

    private void FilterOutliers()
    {
        if(serie.size() <= 2 || !filterOutliers || max-min <= 0 || freqsDict.values().size() < minFreqForOutlierFiltering)
        {
            maxValueNotFiltered = max;
            this.filteredSerie = serie;
            //this.serieDrv1 = drv.GetDerivative(regres, 1);
            //this.serieDrv2 = drv.GetDerivative(serieDrv1, 1);
            return;
        }

        outlierIndexes = new ArrayList<Integer>();
        double meanBias = stat.Mean(bias);
        double stDevBias = stat.StDev(bias);
        if(stDevBias == 0)
        {
            return;
        }
        double maxValue, p;
        int maxIndex = -1;
        boolean[] usedIndexes = new boolean[serie.size()];
        minOutlierFiltered = Double.MAX_VALUE;
        maxValueNotFiltered = -Double.MAX_VALUE;

        while(outlierIndexes.size() < serie.size())
        {
            maxValue = -Double.MAX_VALUE;
            for(int i=0;i<serie.size();i++)
            {
                if(!usedIndexes[i] && serie.get(i).compareTo(maxValue) > 0)
                {
                    maxValue = serie.get(i);
                    maxIndex = i;
                }
            }
            usedIndexes[maxIndex] = true;

            p = ns.pNorm(bias.get(maxIndex), meanBias, stDevBias);

            if(serie.get(maxIndex).compareTo(minNonZeroFreq) > 0 && p > outlierThreshold)
            {
                if(maxValue < minOutlierFiltered)
                {
                    minOutlierFiltered = maxValue;
                }
                outlierIndexes.add(maxIndex);
            }
            else
            {
                maxValueNotFiltered = maxValue;
                Collections.sort(outlierIndexes);
                break;
            }
        }

        if(outlierIndexes.size() > 0)
        {
            results.add(ResultType.OutliersFiltered);
            for(int i=0;i<serie.size();i++)
            {
                if(regres.get(i).compareTo(maxValueNotFiltered) > 0)
                {
                    regres.set(i, maxValueNotFiltered);
                    bias.set(i, 0.0);
                }
                if(serie.get(i).compareTo(maxValueNotFiltered) > 0)
                {
                    serie.set(i, regres.get(i));
                    bias.set(i, 0.0);
                }
            }
            max = maxValueNotFiltered;
        }
        else
        {
            maxValueNotFiltered = max;
        }
        this.filteredSerie = serie;

        //this.serieDrv1 = drv.GetDerivative(regres, 1);
        //this.serieDrv2 = drv.GetDerivative(serieDrv1, 1);
    }

    //endregion

    //region Raw filtering

    private void FilterDailyOutliers(ArrayList<Double> origSerie, int firstIndex)
    {
        int minDailyCount = (int)(1 / (1 - dailyOutlierThreshold));
        double maxFc = 0.5;
        if (origSerie.size() > minDailyCount && fc > (1 - dailyOutlierThreshold))
        {
            if (fc > maxFc)
            {
                this.origSerie = FilterDailyOutliersAll(origSerie, firstIndex);
            }
            else
            {
                this.origSerie = FilterDailyOutliersZerosAndNonZeros(origSerie, firstIndex);
            }
        }
        else
        {
            if(firstIndex > 0)
            {
                this.origSerie = new ArrayList<Double>();
                for (int i = firstIndex; i < origSerie.size();i++)
                {
                    this.origSerie.add(origSerie.get(i));
                }
            }
            else
            {
                this.origSerie = origSerie;
            }
        }
    }

    private ArrayList<Double> FilterDailyOutliersAll(ArrayList<Double> dailySerie, int firstIndex)
    {
        if(dailySerie.size() < minCount)
        {
            return dailySerie;
        }
        double alpha = 1 - this.dailyOutlierThreshold;
        double mean = stat.TrimMeanNonZero(dailySerie, alpha);
        OutlierDetection outlierDetection = new OutlierDetection();
        ArrayList<Double> data = new ArrayList<Double>();
        ArrayList<Double> model = new ArrayList<Double>();
        for(int i=firstIndex;i<dailySerie.size();i++)
        {
            data.add(dailySerie.get(i));
            model.add(mean);
        }
        double lim = outlierDetection.StudDelResLimit(dailySerie, model, this.dailyOutlierThreshold);
        ArrayList<Double> filteredSerie = new ArrayList<Double>();
        //for(int i=0;i<firstIndex;i++) { filteredSerie.Add(dailySerie[i]); }
        ArrayList<Double> studDelRes = outlierDetection.StudDelRes(data, model);
        for(int i=firstIndex;i<dailySerie.size();i++)
        {
            if(studDelRes.get(i-firstIndex).compareTo(lim) < 0 || dailySerie.get(i).compareTo(mean) <= 0)
            {
                filteredSerie.add(dailySerie.get(i));
            }
        }
        return filteredSerie;
    }

    private ArrayList<Double> FilterDailyOutliersZerosAndNonZeros(ArrayList<Double> dailySerie, int firstIndex)
    {
        if(dailySerie.size() < minCount)
        {
            return dailySerie;
        }
        double alpha = 1 - this.dailyOutlierThreshold;
        double mean = stat.TrimMeanNonZero(dailySerie, alpha);
        OutlierDetection outlierDetection = new OutlierDetection();
        ArrayList<Double> data = new ArrayList<Double>();
        ArrayList<Double> dataNonZeros = new ArrayList<Double>();
        ArrayList<Double> model = new ArrayList<Double>();
        ArrayList<Double> modelNonZeros = new ArrayList<Double>();
        for(int i=firstIndex;i<dailySerie.size();i++)
        {
            data.add(dailySerie.get(i));
            if(dailySerie.get(i).compareTo(0) > 0)
            {
                dataNonZeros.add(dailySerie.get(i));
                modelNonZeros.add(mean);
            }
            model.add(mean);
        }
        double lim = outlierDetection.StudDelResLimit(dailySerie, model, this.dailyOutlierThreshold);
        double limNonZeros = outlierDetection.StudDelResLimit(dataNonZeros, modelNonZeros, this.dailyOutlierThreshold);
        ArrayList<Double> filteredSerie = new ArrayList<Double>();
        //for(int i=0;i<firstIndex;i++) { filteredSerie.Add(dailySerie[i]); }
        ArrayList<Double> studDelRes = outlierDetection.StudDelRes(data, model);
        for(int i=firstIndex;i<dailySerie.size();i++)
        {
            if(studDelRes.get(i-firstIndex).compareTo(lim) < 0 || dailySerie.get(i).compareTo(mean) <= 0)
            {
                filteredSerie.add(dailySerie.get(i));
            }
            if(studDelRes.get(i-firstIndex).compareTo(lim) >= 0 && studDelRes.get(i-firstIndex).compareTo(limNonZeros) < 0)
            {
                filteredSerie.add(dailySerie.get(i));
            }
        }
        return filteredSerie;
    }

    //endregion

    //endregion

    //region Normalization

    private void NormalizeCollections()
    {
        if(freqsDict.size() == 1 || serie.size() <= 2)
        {
            normRegres = serie;
            return;
        }
        //serie
        double range = max - min;
        normSerie = new ArrayList<Double>();
        if(range == 0)
        {
            normSerie = serie;
        }
        else
        {
            for(int i=0;i<serie.size();i++)
            {
                normSerie.add(((serie.get(i) - min)/range) * maxNorm);
            }
        }

        //regres
        minReg = Double.MAX_VALUE;
        maxReg = -Double.MAX_VALUE;
        for(double value : serie)
        {
            if(value < minReg)
            {
                minReg = value;
            }
            if(value > maxReg)
            {
                maxReg = value;
            }
        }

        double rangeReg = maxReg - minReg;
        if(range == 0)
        {
            normRegres = regres;
        }
        else
        {
            normRegres = new ArrayList<Double>();
            for(int i=0;i<regres.size();i++)
            {
                normRegres.add(Normalize(regres.get(i), minReg, maxReg, maxNorm));
            }
        }
    }

    private double Normalize(double val, double min, double max, double maxNorm)
    {
        return (val - min)/(max-min) * maxNorm;
    }

    private double UnNormalize(double val, double min, double max, double maxNorm)
    {
        return (val/maxNorm) *(max-min) + min;
    }


    //endregion

    //region Stationarity Methods

    //region Stationarity Clustering Test

    private void HomogeneityClusteringTest(int statPeriodIndex)
    {
        if(serie.size() <= minCount)
        {
            HomogeneityClustering.StatCluster cluster = new HomogeneityClustering.StatCluster();
            SetMinMaxAndFreqs(serie);
            return;
        }
        filteredSerie = new ArrayList<Double>();
        filteredSerie.addAll(serie);
        if (statPeriodIndex > minCountForClustering)
        {
            clusters = hc.Clustering(normRegres, 0, statPeriodIndex-1);
            ArrayList<HomogeneityClustering.StatCluster> newClusters = hc.Clustering(normRegres, statPeriodIndex, serie.size() - 1);
            clusters.addAll(newClusters);
        }
        else
        {
            clusters = hc.Clustering(normRegres, 0, serie.size()-1);
        }
        means = hc.GetMeans(clusters);
        stDevs = hc.GetStDevs(clusters);

        int i = clusters.size()-1;
        firstCluster = new HomogeneityClustering.StatCluster(clusters.get(i));

        ArrayList<HomogeneityClustering.StatCluster> walls = new ArrayList<HomogeneityClustering.StatCluster>();

        if(testType == HomogeneityClustering.TestType.Statistics)
        {
            while (serie.size() - firstCluster.FirstIndex < minCountForClustering && i > 0)
            {
                firstCluster.Join(clusters.get(--i));
            }
            while(i>0)
            {

                if(clusters.get(i-1).Mean > firstCluster.Mean && clusters.get(i-1).size() <= minWall)
                {
                    walls.add(clusters.get(i-1));
                    i--;
                    continue;
                }
                else if (clusters.get(i - 1).Mean <= firstCluster.Mean && clusters.get(i - 1).size() <= minCountForClustering)
                {
                    if(firstCluster.Mean - clusters.get(i-1).Mean <= hc.OrigMaxDiff)
                    {
                        firstCluster.Join(clusters.get(i-1));
                    }
                    i--;
                    continue;

                }
                else
                {
                    double maxThresholdForAdvance = hc.OrigMaxDiff;

                    if(Math.abs(clusters.get(i-1).Mean - firstCluster.Mean) < maxThresholdForAdvance * weights[firstCluster.FirstIndex] && (firstCluster.Mean == 0 || clusters.get(i-1).Mean / firstCluster.Mean < 2))
                {
                        firstCluster.Join(clusters.get(i-1));
                        i--;
                    }
                    else
                    {
                        break;
                    }
                }
            }
        }
        else
        {
            while (i > 0 && (serie.size() - firstCluster.FirstIndex < minCountForClustering || hc.AreHomogeneous(clusters.get(i - 1), firstCluster)))
            {
                if(clusters.get(i-1).Mean > firstCluster.Mean && clusters.get(i-1).size() <= 2)
                {
                    walls.add(clusters.get(i-1));
                }
                firstCluster.Join(clusters.get(i-1));
                i--;
            }
        }

        if (firstCluster.FirstIndex > firstIndexToFcst)
        {
            firstIndexToFcst = firstCluster.FirstIndex;

        }
        if (serie.size() - firstIndexToFcst < minCountForClustering)
        {
            firstIndexToFcst = serie.size() - minCount - 1;
        }

        //filtrado walls
        if(outlierIndexes == null)
        {
            outlierIndexes = new ArrayList<Integer>();
        }
        ArrayList<Double> trunkMeans = new ArrayList<Double>();
        for(int j=firstIndexToFcst;j<means.size();j++)
        {
            trunkMeans.add(means.get(j));
        }
        double meanTrunk = stat.Mean(trunkMeans);
        double stDevTrunk = stat.StDev(trunkMeans);
        if(stDevTrunk == 0)
        {
            return;
        }
        double p;
        boolean clustFiltered = false;
        for(int j=0;j<walls.size();j++)
        {
            p = ns.pNorm(walls.get(j).Mean, meanTrunk, stDevTrunk);
            if(p > outlierThreshold)
            {
                if(walls.get(j).Mean < minOutlierFiltered)
                {
                    minOutlierFiltered = walls.get(j).Mean;
                    clustFiltered = true;
                }
                for(int k=walls.get(j).FirstIndex;k<=walls.get(j).LastIndex;k++)
                {
                    outlierIndexes.add(k);
                    double fixedValue = UnNormalize(firstCluster.Mean, minReg, maxReg, maxNorm);
                    serie.set(k, fixedValue);
                    regres.set(k, fixedValue);
                    bias.set(k, 0);
                    means.set(k, fixedValue);
                }
            }
        }
        if(outlierIndexes.size() > 1)
        {
            Collections.sort(outlierIndexes);
        }
        if(clustFiltered)
        {
            results.add(ResultType.ClusteringFiltered);
        }

        //serie.RemoveRange(0, firstIndexToFcst);
        //regres.RemoveRange(0, firstIndexToFcst);
        //bias.RemoveRange(0, firstIndexToFcst);
        SetMinMaxAndFreqs(serie);
        statPeriodMean = UnNormalize(firstCluster.Mean, min, max, maxNorm);
    }

    private void SetMinMaxAndFreqs(ArrayList<Double> serie)
    {
        min = Double.MAX_VALUE;
        max = -Double.MAX_VALUE;
        freqsDict.clear();
        for(double value : serie)
        {
            if(value < min)
            {
                min = Math.ceil(value);
            }
            if(value > max)
            {
                max = Math.ceil(value);
            }
            AddFreq(value);
        }
    }

    //endregion 

    //region Weighting Clustering 

    private void SetClusterWeights()
    {
        if(serie.size() < minCount || clusters == null || clusters.isEmpty())
        {
            return;
        }
        double minPVal = Double.MAX_VALUE;
        int minPValIndex = -1;
        int minPValIndex2 = -1;
        double pVal;
        for(int i=0;i<clusters.size();i++)
        {
            if(pValue > 0)
            {
                pVal = hypoTest.MannWhitney(serie, clusters.get(i).FirstIndex, true);
                if(pVal < minPVal)
                {
                    minPVal = pVal;
                    minPValIndex = i;
                }
                if(pVal == minPVal)
                {
                    minPValIndex2 = i;
                }
            }
            if(clusters.get(i).FirstIndex >= firstIndexToFcst)
            {
                clusters.get(i).Weight = 1.0;
            }
            else
            {
                clusters.get(i).Weight = CalcWeigth(clusters.get(i));
            }

            for (int j = clusters.get(i).FirstIndex; j <= clusters.get(i).LastIndex; j++)
            {
                weights[j] *= clusters.get(i).Weight;
            }
        }
    }

    private double CalcWeigth(HomogeneityClustering.StatCluster sc)
    {
        //normal intersection : 
        //double w1 = firstCluster.Count;
        //double w2 = sc.Count;
        double w1 = ((double)firstCluster.size()/(double)(firstCluster.size() + sc.size())) * maxNorm;
        double w2 = ((double)sc.size() / (double)(firstCluster.size() + sc.size())) * maxNorm;
        double m1 = firstCluster.Mean;
        double m2 = sc.Mean;
        double s1 = firstCluster.StDev;
        double s2 = sc.StDev;

        double weight = ns.GetNormalIntersection(m1, s1, w1, m2, s2, w2, true);
        //double weight = hypoTest.MannWhitney(firstCluster.Values, sc.Values, true);

        if (rootExpForClusterWeighting > 1)
        {
            weight = Math.pow(weight, (1.0 / rootExpForClusterWeighting));
        }
        return weight;
    }

    //endregion

    //endregion

    //region Weighting Methods

    private void CalculateWeights()
    {
        if (lag < nSpans * 2)
        {
            weights = func.CalcTempWeights(serie.size(), forgetInitPeriods, forgetEndPeriods, forgetEndProportion);
        }
        else
        {
            weights = CalculateSeasonalWeights(serie, lag);
        }
    }

    private double[] CalculateSeasonalWeights(ArrayList<Double> serie, int lag)
    {
        double[] wghV = func.CalcTempWeights(lag-nSpans, 0, (lag-nSpans)/2, forgetEndProportion);
        for (int i = 0; i < lag / 2; i++)
        {
            wghV[i] = wghV[wghV.length - i-1];
        }

        double[] wgh1 = new double[nSpans];
        for (int i = 0; i < nSpans; i++)
        {
            wgh1[i] = 1.00;
        }

        ArrayList<Double> weightList = new ArrayList<Double>();
        while (weightList.size() < serie.size())
        {
            tangible.DoubleLists.addPrimitiveArrayToList(wgh1, weightList);
            tangible.DoubleLists.addPrimitiveArrayToList(wghV, weightList);
        }
        while (weightList.size() > serie.size())
        {
            weightList.remove(0);
        }
        weights = tangible.DoubleLists.toArray(weightList);
        return weights;
    }

    //endregion
    
    import java.util.*;

    public class SignalCharacterization2
    {

    		//region Auxiliar Methods

    		private double GetLoadFactor(ArrayList<Double> serie)
    		{
    			int tot = 0;
    			for (double val : serie)
    			{
    				if (val > 0)
    				{
    					tot++;
    				}
    			}

    			double lfactor = (double)tot / ((double)serie.size());
    			return lfactor;
    		}

    		private double GetWeightedLoadFactor(ArrayList<Double> serie)
    		{
    			int forgetInitIndex = serie.size() - Math.min(serie.size(), forgetInitPeriods);
    			int forgetEndIndex = serie.size() - Math.min(serie.size(), forgetEndPeriods);

    			double n = (double)serie.size();

    			double w1 = (serie.size() - forgetInitIndex) * 1;
    			double w2 = (forgetInitIndex - forgetEndIndex) * ((1 - forgetEndProportion) / 2);
    			double w3 = forgetEndIndex * forgetEndProportion;

    			double totW = w1 + w2 + w3;
    			w1 /= totW;
    			w2 /= totW;
    			w3 /= totW;
    			double lf1 = 0;
    			double lf2 = 0;
    			double lf3 = 0;
    			double tot1 = 0;
    			for (int i = forgetInitIndex;i < serie.size();i++)
    			{
    				if (serie.get(i).compareTo(0) > 0)
    				{
    					tot1++;
    				}
    			}
    			if (serie.size() - forgetInitIndex > 0)
    			{
    				lf1 = tot1 / (double)(serie.size() - forgetInitIndex);
    			}

    			double tot2 = 0;
    			for (int i = forgetEndIndex; i < forgetInitIndex; i++)
    			{
    				if (serie.get(i).compareTo(0) > 0)
    				{
    					tot2++;
    				}
    			}
    			if (forgetInitIndex - forgetEndIndex > 0)
    			{
    				lf2 = tot2 / (double)(forgetInitIndex - forgetEndIndex);
    			}

    			double tot3 = 0;
    			for (int i = 0; i < forgetEndIndex; i++)
    			{
    				if (serie.get(i).compareTo(0) > 0)
    				{
    					tot3++;
    				}
    			}
    			if (forgetEndIndex > 0)
    			{
    				lf3 = tot3 / (double)(forgetEndIndex);
    			}


    			 double lfactor = lf1 * w1 + lf2 * w2 + lf3 * w3;
    			 return lfactor;
    		}

    		private ArrayList<Double> GetResiduals(ArrayList<Double> serie)
    		{
    			Result res = poly.LSRegression(serie);
    			ArrayList<Double> coeffs = new ArrayList<Double>();
    			coeffs.add(res.c0);
    			coeffs.add(res.c1);
    			ArrayList<Double> x = new ArrayList<Double>();
    			for (int i = 0;i < serie.size();i++)
    			{
    				x.add((double)i);
    			}
    			ArrayList<Double> lsReg = poly.GetRegValues(x, coeffs);
    			ArrayList<Double> resids = func.Substract(serie, lsReg, true, false);
    			return resids;
    		}

    		private void AddFreq(double val)
    		{
    			if (val > 0 && val < minNonZeroFreq)
    			{
    				minNonZeroFreq = val;
    			}
    			int key = (int)Math.ceil(val);
    			if (!freqsDict.ContainsKey(key))
    			{
    				freqsDict.Add(key, 1.0);
    			}
    			else
    			{
    				freqsDict[key] = freqsDict[key] + 1.0;
    			}
    		}

    		private void SubstractFreq(double val)
    		{
    			int key = (int)Math.ceil(val);
    			if (freqsDict.ContainsKey(key))
    			{
    				freqsDict[key] = freqsDict[key] - 1.0;
    			}
    			if (freqsDict[key] == 0)
    			{
    				freqsDict.Remove(key);
    			}
    		}

    		private ArrayList<Integer> GetFiniteJumpCandidates(ArrayList<Double> serie)
    		{
    			return poly.ZerosNegPosIndexes(serieDrv2);
    		}

    		private void AddInitialValues(ArrayList<Double> serie, ArrayList<Double> initialValues, int ini, int fin)
    		{
    			for (int i = fin;i >= ini;i--)
    			{
    				serie.add(0, initialValues.get(i));
    			}
    		}

    		//endregion

    		//endregion

    		//region Enums

    		/**  Result of calculation 
    		*/
    		public enum ResultType
    		{

    			/**  Some outliers have been filtered (at least one) 
    			*/
    			OutliersFiltered,

    			/**  Total periods is below the minimum allowed 
    			*/
    			LessThanMinPeriods,

    			/**  The time series is stationary 
    			*/
    			Stationary,

    			/**  Some periods have been filtered by clustering processes 
    			*/
    			ClusteringFiltered,

    			/**  Sparse and only one cluster as a result of clustering process 
    			*/
    			Sparse,

    			/**  Treated as sparse 
    			*/
    			AsSparse;

    			public int getValue()
    			{
    				return this.ordinal();
    			}

    			public static ResultType forValue(int value)
    			{
    				return values()[value];
    			}
    		}

    		//endregion

    }
}

