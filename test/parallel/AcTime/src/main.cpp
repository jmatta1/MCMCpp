#include<iostream>
//#include<fstream>
//#include<iomanip>
#include"Movers/Diagnostic/AutoRegressiveMove.h"
#include"Analysis/AutoCorrCalc.h"
#include"ParallelEnsembleSampler.h"

using MCMC::ParallelEnsembleSampler;
namespace Mover=MCMC::Mover;
namespace Analysis=MCMC::Analysis;

//std::string FileNames[5] = {"p1.csv", "p2.csv", "p3.csv", "p4.csv", "p5.csv"};

int main()
{
    const int runNumber = 0;
    const int numWalkers = 240;
    const int numParams = 5;
    const int numSteps = 131071;
    const int numThreads = 3;
    double offsets[numParams] = {1.0, 1.0, 1.0, 1.0, 1.0};
    double phis[numParams] = {0.8,  // AC = 9
                              0.904761904762,  // AC = 20
                              0.9354838709677,  // AC = 30
                              0.9672131147541,  // AC = 60
                              0.990050200903734685};  // AC = 200
    double vars[numParams] = {1.0, 1.0, 1.0, 1.0, 1.0};
    
    std::cout<<"Building initial mover"<<std::endl;
    Mover::AutoRegressiveMove<double> mover(numParams, offsets, phis, vars);
    
    std::cout<<"Building parallel ensemble sampler"<<std::endl;
    ParallelEnsembleSampler<double, Mover::AutoRegressiveMove<double> > sampler(runNumber, numThreads, numWalkers, numParams, mover);
    
    std::cout<<"Getting initial values"<<std::endl;
    double* initVals = new double[numWalkers*numParams];
    double* auxVals = new double[numWalkers];
    mover.getInitialPoints(initVals, auxVals, numWalkers);
    
    std::cout<<"Setting initial values"<<std::endl;
    sampler.setInitialWalkerPos(initVals, auxVals);
    delete[] initVals;
    delete[] auxVals;
    
    std::cout<<"Running parallel ensemble sampler for autoregressive model"<<std::endl;
    sampler.runMCMC(numSteps);
    
    /*std::cout<<"Writing out chains"<<std::endl;
    auto end = sampler.getStepIttEnd();
    for(int j=0; j<numParams; ++j)
    {
        std::ofstream output(FileNames[j]);
        for(auto itt = sampler.getStepIttBegin(); itt != end; ++itt)
        {
            
            output<<(*itt)[j];
            for(int i=1; i<numWalkers; ++i)
            {
                output<<", "<<std::scientific<<std::setprecision(16)<<(*itt)[i*numParams+j];
            }
            output<<'\n';
        }
        output<<std::flush;
    }*/
    
    std::cout<<"Calculating integrated autocorrelation times"<<std::endl;
    Analysis::AutoCorrCalc<double> acCalc(numParams, numWalkers);
    auto startItt = sampler.getStepIttBegin();
    auto endItt = sampler.getStepIttEnd();
    acCalc.calcAutoCorrTimes(startItt, endItt, sampler.getStoredSteps());
    
    std::cout<<"P0 Calculated AutoCorrelation Time: "<<acCalc.retrieveAutoCorrelationTime(0)<<" - Actual:   9"<<std::endl;
    std::cout<<"P1 Calculated AutoCorrelation Time: "<<acCalc.retrieveAutoCorrelationTime(1)<<" - Actual:  20"<<std::endl;
    std::cout<<"P2 Calculated AutoCorrelation Time: "<<acCalc.retrieveAutoCorrelationTime(2)<<" - Actual:  30"<<std::endl;
    std::cout<<"P3 Calculated AutoCorrelation Time: "<<acCalc.retrieveAutoCorrelationTime(3)<<" - Actual:  60"<<std::endl;
    std::cout<<"P4 Calculated AutoCorrelation Time: "<<acCalc.retrieveAutoCorrelationTime(4)<<" - Actual: 200"<<std::endl;
}
