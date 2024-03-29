#include<iostream>
#include<random>
#include<cmath>
#include<string>
#include<sstream>
#include<fstream>
#include"Common/SkewedGaussian.h"
#include"Utility/pcg-cpp/include/pcg_random.hpp"
#include"Analysis/AutoCorrCalc.h"
#include"Analysis/CovarianceMatrix.h"
#include"Analysis/CornerHistograms.h"
#include"Analysis/PercentileAndMaximumFinder.h"
#include"Movers/MetropolisHastings.h"
#include"EnsembleSampler.h"
using MCMC::EnsembleSampler;
namespace Mover=MCMC::Mover;
namespace Analysis=MCMC::Analysis;

const char* HelpStr=" <0|1|2> [no_write]\n  Where 0 means use the ideal "
              "covariance matrix\n        1 means use the identity covariance "
              "matrix\n    and 2 means use the bad covariance matrix\n        "
              "Using [no_write] prevents writing the chains";

template<class ParamType>
void generateInitialValues(ParamType* initVals, ParamType* auxVals, SkewedGaussianTwoDim<ParamType>& likelihood,
                           int numWalkers, int numParams, int extraRunNumber);

int main(int argc, char* argv[])
{
    typedef Mover::MetropolisHastings<double, SkewedGaussianTwoDim<double> > Walker;
    int option=0;
    bool writeChains = true;

    // check to make sure that the user entered the correct number of arguments
    if(argc != 2 && argc != 3)
    {
        std::cout<<"Usage:\n    "<<argv[0]<<HelpStr<<std::endl;
        return 1;
    }

    std::string argOne(argv[1]);

    // check to make sure that the second argument is 0, 1, or 2
    if(argOne != "0" && argOne != "1" && argOne != "2")
    {
        std::cout<<"Bad first argument: \""<<argv[1]<<"\""<<std::endl;
        std::cout<<"Usage:\n    "<<argv[0]<<HelpStr<<std::endl;
        return 1;
    }
    std::istringstream converter;
    converter.str(argv[1]);
    converter >> option;

    // if we got three arguments then check if the third one is correct
    if(argc == 3)
    {
        std::string argTwo(argv[2]);
        if(argTwo != "no_write")
        {
            std::cout<<"Bad second argument: \""<<argv[2]<<"\""<<std::endl;
            std::cout<<"Usage:\n    "<<argv[0]<<HelpStr<<std::endl;
            return 1;
        }
        writeChains = false;
    }

    const int runNumber = 0;
    const int extraRunNumber = 53;
    const int numWalkers = 320;
    const int numParams = 2;
    const int numSteps = 40019;
    const int cornerBinning = 100;
    const double eps = 0.13;
    const double idealCovar[4] = {(1.0+eps), (1.0-eps)/2.0, (1.0-eps)/2.0, (1.0+eps)/4.0};
    const double badCovar[4]   = {(2.0/3.0), -5.0, -5.0, 50.0};
    //const double identCovar[4] = {1.0, 0.0, 0.0, 1.0};
    
    std::cout<<"Building Custom Distribution"<<std::endl;
    SkewedGaussianTwoDim<double> likelihood(eps);
    
    std::cout<<"Building initial mover"<<std::endl;
    Walker mover(numParams, runNumber, likelihood);
    if(option==0)
    {
        std::cout<<"Setting Metropolis Hastings sample covariance matrix to ideal case"<<std::endl;
        mover.setCovarMat(idealCovar);
    }
    else if(option==1)
    {
        std::cout<<"Setting Metropolis Hastings sample covariance matrix to identity case"<<std::endl;
        mover.setIdentityCovarMat();
    }
    else if(option==2)
    {
        std::cout<<"Setting Metropolis Hastings sample covariance matrix to bad case"<<std::endl;
        mover.setCovarMat(badCovar);
    }
    
    std::cout<<"Building sampler"<<std::endl;
    EnsembleSampler<double, Walker> sampler(runNumber, numWalkers, numParams, mover);
    
    std::cout<<"Setting sampler to skip 9 points for every point it remembers"<<std::endl;
    sampler.setSlicingMode(true, 10);
    
    std::cout<<"Generating initial values"<<std::endl;
    double* initVals = new double[numWalkers*numParams];
    double* auxVals = new double[numWalkers];
    generateInitialValues(initVals, auxVals, likelihood, numWalkers, numParams, extraRunNumber);
    
    std::cout<<"Setting initial values"<<std::endl;
    sampler.setInitialWalkerPos(initVals, auxVals);
    std::cout<<"Deleting initial value arrays"<<std::endl;
    delete[] initVals;
    delete[] auxVals;
    
    std::cout<<"Running ensemble sampler"<<std::endl;
    bool status = sampler.runMCMC(numSteps);
    
    if(status)
    {
        std::cout<<"Sampling completed normally"<<std::endl;
    }
    else
    {
        std::cout<<"Sampling finished when chain ran out of space"<<std::endl;
    }
    std::cout<<"Discarding 20 points for burn in."<<std::endl;
    sampler.sliceAndBurnChain(1, 20);
    std::cout<<"\nAcceptance Fraction: "<<sampler.getAcceptedSteps()<<"/"
             <<sampler.getTotalSteps()<<" | "<<sampler.getAcceptanceFraction()
             <<std::endl;
    if(writeChains)
    {
        std::cout<<"\nWriting out chains"<<std::endl;
        std::ofstream output;
        if(option==0)
        {
            output.open("linearChainsIdeal.csv");
        }
        else if(option==1)
        {
            output.open("linearChainsIdentity.csv");
        }
        else if(option==2)
        {
            output.open("linearChainsBad.csv");
        }
        auto end = sampler.getParamSetIttEnd();
        for(auto itt = sampler.getParamSetIttBegin(); itt != end; ++itt)
        {
            output<<(*itt)[0]<<", "<<(*itt)[1]<<std::endl;
        }
        output.close();
    }
    std::cout<<"\nCalculating integrated autocorrelation times"<<std::endl;
    Analysis::AutoCorrCalc<double> acCalc(numParams, numWalkers);
    auto startItt = sampler.getStepIttBegin();
    auto endItt = sampler.getStepIttEnd();
    acCalc.calcAutoCorrTimes(startItt, endItt, sampler.getStoredSteps());
    
    double p0Ac = acCalc.retrieveAutoCorrelationTime(0);
    double p1Ac = acCalc.retrieveAutoCorrelationTime(1);
    
    std::cout<<"P0 Calculated AutoCorrelation Time: "<<p0Ac<<std::endl;
    std::cout<<"P1 Calculated AutoCorrelation Time: "<<p1Ac<<std::endl;
    
    std::cout<<"\nCalculating the covariance matrix with slicing"<<std::endl;
    Analysis::CovarianceMatrix<double> cmCalc(numParams, numWalkers);
    int sliceInterval = static_cast<int>((p0Ac<p1Ac)?std::ceil(p0Ac):std::ceil(p1Ac));
    cmCalc.calculateCovar(startItt, endItt, sliceInterval);
    std::cout<<"Covariance matrix with slicing"<<std::endl;
    std::cout<<cmCalc.getCovarianceMatrixElement(0, 0)<<", "<<cmCalc.getCovarianceMatrixElement(0, 1)<<"\n";
    std::cout<<cmCalc.getCovarianceMatrixElement(1, 0)<<", "<<cmCalc.getCovarianceMatrixElement(1, 1)<<"\n";
    std::cout<<"Correlation matrix with slicing"<<std::endl;
    std::cout<<cmCalc.getCorrelationMatrixElement(0, 0)<<", "<<cmCalc.getCorrelationMatrixElement(0, 1)<<"\n";
    std::cout<<cmCalc.getCorrelationMatrixElement(1, 0)<<", "<<cmCalc.getCorrelationMatrixElement(1, 1)<<"\n";
    
    std::cout<<"\nGenerating Corner Histograms"<<std::endl;
    Analysis::CornerHistograms<double> cornerHists(numParams, numWalkers, cornerBinning);
    cornerHists.calculateHistograms(startItt, endItt);
    std::cout<<"Writing Corner Histograms"<<std::endl;
    cornerHists.saveHistsCsvFormat("chainHist");
    
    std::cout<<"\nGenerating Percentile and Peak Finding Histograms"<<std::endl;
    Analysis::PercentileAndMaximumFinder<double> pamf(numParams, numWalkers, 5*cornerBinning);
    pamf.processChainData(startItt, endItt, 1);
    std::cout<<"Writing High-Res Histograms"<<std::endl;
    pamf.writeHistogramsInCsvFormat("percentileHistograms");
    
    std::cout<<"Finding peaks and errors"<<std::endl;
    double peak0 = pamf.getValueOfPeak(0);
    double percentile0 = pamf.getPercentileFromValue(0, peak0);
    double loPeak0 = pamf.getValueFromPercentile(0, percentile0-34.1);
    double hiPeak0 = pamf.getValueFromPercentile(0, percentile0+34.1);
    double peak1 = pamf.getValueOfPeak(1);
    double percentile1 = pamf.getPercentileFromValue(1, peak1);
    double loPeak1 = pamf.getValueFromPercentile(1, percentile1-34.1);
    double hiPeak1 = pamf.getValueFromPercentile(1, percentile1+34.1);
    std::cout<<"The percentiles, values and errors of the peaks are:\n"
             <<"Parameter, Peak Percentile, Peak Value, 34.1% down, 34.1% up"<<std::endl;
    std::cout<<"P0: "<<percentile0<<", "<<peak0<<", "<<loPeak0<<", "<<hiPeak0<<std::endl;
    std::cout<<"P1: "<<percentile1<<", "<<peak1<<", "<<loPeak1<<", "<<hiPeak1<<std::endl;
    
    std::cout<<"Finding percentiles"<<std::endl;
    double cv0 = pamf.getValueFromPercentile(0, 50);
    double lv0 = pamf.getValueFromPercentile(0, 15.9);
    double hv0 = pamf.getValueFromPercentile(0, 84.1);
    double cv1 = pamf.getValueFromPercentile(1, 50);
    double lv1 = pamf.getValueFromPercentile(1, 15.9);
    double hv1 = pamf.getValueFromPercentile(1, 84.1);
    std::cout<<"The 15.9, 50, 84.1 percentiles are:"<<std::endl;
    std::cout<<"P0: "<<lv0<<", "<<cv0<<", "<<hv0<<std::endl;
    std::cout<<"P1: "<<lv1<<", "<<cv1<<", "<<hv1<<std::endl;
    
    std::cout<<"\nShutting down"<<std::endl;
    
    return 0;
}

template<class ParamType>
void generateInitialValues(ParamType* initVals, ParamType* auxVals, SkewedGaussianTwoDim<ParamType>& likelihood, int numWalkers, int numParams, int extraRunNumber)
{
    const double averages[2] = {3.1, -3.5};
    const double deviations[2] = {3.3, 4.1};
    pcg32 engine(extraRunNumber);
    std::normal_distribution<double> nd1(averages[0], deviations[0]);
    std::normal_distribution<double> nd2(averages[1], deviations[1]);
    std::normal_distribution<double> nd3(averages[0], deviations[1]);
    std::normal_distribution<double> nd4(averages[1], deviations[0]);
    for(int i=0; i<numWalkers; i+=16)
    {
        initVals[numParams*(i+0)] = nd1(engine);
        initVals[numParams*(i+0)+1] = nd1(engine);
        auxVals[(i+0)]=likelihood.calcLogPostProb(initVals + numParams*(i+0));
        initVals[numParams*(i+1)] = nd1(engine);
        initVals[numParams*(i+1)+1] = nd2(engine);
        auxVals[(i+1)]=likelihood.calcLogPostProb(initVals + numParams*(i+1));
        initVals[numParams*(i+2)] = nd1(engine);
        initVals[numParams*(i+2)+1] = nd3(engine);
        auxVals[(i+2)]=likelihood.calcLogPostProb(initVals + numParams*(i+2));
        initVals[numParams*(i+3)] = nd1(engine);
        initVals[numParams*(i+3)+1] = nd4(engine);
        auxVals[(i+3)]=likelihood.calcLogPostProb(initVals + numParams*(i+3));
        
        initVals[numParams*(i+4)] = nd2(engine);
        initVals[numParams*(i+4)+1] = nd1(engine);
        auxVals[(i+4)]=likelihood.calcLogPostProb(initVals + numParams*(i+4));
        initVals[numParams*(i+5)] = nd2(engine);
        initVals[numParams*(i+5)+1] = nd2(engine);
        auxVals[(i+5)]=likelihood.calcLogPostProb(initVals + numParams*(i+5));
        initVals[numParams*(i+6)] = nd2(engine);
        initVals[numParams*(i+6)+1] = nd3(engine);
        auxVals[(i+6)]=likelihood.calcLogPostProb(initVals + numParams*(i+6));
        initVals[numParams*(i+7)] = nd2(engine);
        initVals[numParams*(i+7)+1] = nd4(engine);
        auxVals[(i+7)]=likelihood.calcLogPostProb(initVals + numParams*(i+7));
        
        initVals[numParams*(i+8)] = nd3(engine);
        initVals[numParams*(i+8)+1] = nd1(engine);
        auxVals[(i+8)]=likelihood.calcLogPostProb(initVals + numParams*(i+8));
        initVals[numParams*(i+9)] = nd3(engine);
        initVals[numParams*(i+9)+1] = nd2(engine);
        auxVals[(i+9)]=likelihood.calcLogPostProb(initVals + numParams*(i+9));
        initVals[numParams*(i+10)] = nd3(engine);
        initVals[numParams*(i+10)+1] = nd3(engine);
        auxVals[(i+10)]=likelihood.calcLogPostProb(initVals + numParams*(i+10));
        initVals[numParams*(i+11)] = nd3(engine);
        initVals[numParams*(i+11)+1] = nd4(engine);
        auxVals[(i+11)]=likelihood.calcLogPostProb(initVals + numParams*(i+11));
        
        initVals[numParams*(i+12)] = nd4(engine);
        initVals[numParams*(i+12)+1] = nd1(engine);
        auxVals[(i+12)]=likelihood.calcLogPostProb(initVals + numParams*(i+12));
        initVals[numParams*(i+13)] = nd4(engine);
        initVals[numParams*(i+13)+1] = nd2(engine);
        auxVals[(i+13)]=likelihood.calcLogPostProb(initVals + numParams*(i+13));
        initVals[numParams*(i+14)] = nd4(engine);
        initVals[numParams*(i+14)+1] = nd3(engine);
        auxVals[(i+14)]=likelihood.calcLogPostProb(initVals + numParams*(i+14));
        initVals[numParams*(i+15)] = nd4(engine);
        initVals[numParams*(i+15)+1] = nd4(engine);
        auxVals[(i+15)]=likelihood.calcLogPostProb(initVals + numParams*(i+15));
    }
}
