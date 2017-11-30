#include<iostream>
#include"Movers/Diagnostic/AutoRegressiveMove.h"
#include"Analysis/AutoCorrCalc.h"
#include"EnsembleSampler.h"

using MCMC::EnsembleSampler;
namespace Mover=MCMC::Mover;
namespace Analysis=MCMC::Analysis;

int main()
{
    const int runNumber = 0;
    const int numWalkers = 100;
    const int numParams = 5;
    const int numSteps = 262143;
    double offsets[numParams] = {1.0, 1.0, 1.0, 1.0, 1.0};
    double phis[numParams] = {0.8,  // AC = 9
                             0.904761904762,  // AC = 20
                             0.9354838709677,  // AC = 30
                             0.9672131147541,  // AC = 60
                             0.990050200903734685};  // AC = 200
    double vars[numParams] = {1.0, 1.0, 1.0, 1.0, 1.0};
    
    std::cout<<"Building initial mover"<<std::endl;
    Mover::AutoRegressiveMove<double> mover(numParams, offsets, phis, vars);
    
    std::cout<<"Getting initial values"<<std::endl;
    double* initVals = new double[numWalkers*numParams];
    double* auxVals = new double[numWalkers];
    mover.getInitialPoints(initVals, auxVals, numWalkers);
    
    std::cout<<"Building sampler"<<std::endl;
    EnsembleSampler<double, Mover::AutoRegressiveMove<double> > sampler(runNumber, numWalkers, numParams, mover);
    std::cout<<"Setting initial values"<<std::endl;
    sampler.setInitialWalkerPos(initVals, auxVals);
    delete[] initVals;
    delete[] auxVals;
    
    std::cout<<"Running ensemble sampler for autoregressive model"<<std::endl;
    sampler.runMCMC(numSteps);
    
    std::cout<<"Calculating integrated autocorrelation times"<<std::endl;
    Analysis::AutoCorrCalc<double> acCalc(numParams, numWalkers);
    auto startItt = sampler.getStepIttBegin();
    auto endItt = sampler.getStepIttEnd();
    
    //num steps +1 to include the init location
    acCalc.calcAutoCorrTimes(startItt, endItt, numSteps+1);
    std::cout<<"P0 Calculated AutoCorrelation Time: "<<acCalc.retrieveAutoCorrelationTime(0)<<" Actual:   9"<<std::endl;
    std::cout<<"P1 Calculated AutoCorrelation Time: "<<acCalc.retrieveAutoCorrelationTime(1)<<" Actual:  20"<<std::endl;
    std::cout<<"P2 Calculated AutoCorrelation Time: "<<acCalc.retrieveAutoCorrelationTime(2)<<" Actual:  30"<<std::endl;
    std::cout<<"P3 Calculated AutoCorrelation Time: "<<acCalc.retrieveAutoCorrelationTime(3)<<" Actual:  60"<<std::endl;
    std::cout<<"P4 Calculated AutoCorrelation Time: "<<acCalc.retrieveAutoCorrelationTime(4)<<" Actual: 200"<<std::endl;
}
