#include<iostream>
#include"Movers/Diagnostic/AutoRegressiveMove.h"
#include"Analysis/AutoCorrCalc.h"
#include"ParallelEnsembleSampler.h"

using MCMC::ParallelEnsembleSampler;
namespace Mover=MCMC::Mover;
namespace Analysis=MCMC::Analysis;

int main()
{
    //when setting the run number to zero (and thus the random seed to zero, and using 4 threads) the
    //following values are obtained as auto correlation times for each chain
    //phi = 0.8 gives a calculated value of: 9.01951
    //phi = 0.904761904762 gives a calculated value of: 19.9437
    //phi = 0.9354838709677 gives a calculated value of: 29.7831
    //phi = 0.9672131147541 gives a calculated value of: 59.8488
    //phi = 0.990050200903734685 gives a calculated value of: 196.85 (I suspect truncation of the value when translated into a literal)
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
    
    std::cout<<"Calculating integrated autocorrelation times"<<std::endl;
    Analysis::AutoCorrCalc<double> acCalc(numParams, numWalkers);
    auto startItt = sampler.getStepIttBegin();
    auto endItt = sampler.getStepIttEnd();
    acCalc.calcAutoCorrTimes(startItt, endItt, sampler.getStoredSteps());
    
    std::cout<<"P0 Calculated AutoCorrelation Time:   "<<acCalc.retrieveAutoCorrelationTime(0)<<" Expected:   9.01951 - Actual:   9"<<std::endl;
    std::cout<<"P1 Calculated AutoCorrelation Time:  "<<acCalc.retrieveAutoCorrelationTime(1)<<"  Expected:  19.9437  - Actual:  20"<<std::endl;
    std::cout<<"P2 Calculated AutoCorrelation Time:  "<<acCalc.retrieveAutoCorrelationTime(2)<<"  Expected:  29.7831  - Actual:  30"<<std::endl;
    std::cout<<"P3 Calculated AutoCorrelation Time:  "<<acCalc.retrieveAutoCorrelationTime(3)<<"  Expected:  59.8488  - Actual:  60"<<std::endl;
    std::cout<<"P4 Calculated AutoCorrelation Time: "<<acCalc.retrieveAutoCorrelationTime(4)<<"    Expected: 196.85    - Actual: 200"<<std::endl;
}
