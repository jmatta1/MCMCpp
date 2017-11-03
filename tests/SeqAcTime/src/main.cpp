#include<iostream>
#include"Movers/AutoRegressiveMover.h"
#include"Analysis/AutoCorrCalc.h"
#include"EnsembleSampler.h"
using MarkovChainMonteCarlo::EnsembleSampler;
namespace Mover=MarkovChainMonteCarlo::Mover;
namespace Analysis=MarkovChainMonteCarlo::Analysis;

int main()
{
    const int numWalkers = 100;
    const int numParams = 5;
    const int numSteps = 200000;
    float offsets[numParams] = {1.0f, 1.0f, 1.0f, 1.0f, 1.0f};
    float phis[numParams] = {0.8f,  // AC = 9
                             0.904761904762f,  // AC = 20
                             0.9354838709677f,  // AC = 30
                             0.9672131147541f,  // AC = 60
                             0.990050200903734685f};  // AC = 200
    float vars[numParams] = {1.0f, 1.0f, 1.0f, 1.0f, 1.0f};
    
    std::cout<<"Building initial mover"<<std::endl;
    Mover::AutoRegressiveMove<float> mover(numParams, offsets, phis, vars);
    
    std::cout<<"Getting initial values"<<std::endl;
    float* initVals = new float[numWalkers*numParams];
    float* auxVals = new float[numWalkers];
    mover.getInitialPoints(initVals, auxVals, numWalkers);
    
    std::cout<<"Building sampler"<<std::endl;
    EnsembleSampler<float, Mover::AutoRegressiveMove<float> > sampler(0, numWalkers, numParams, mover);
    std::cout<<"Setting initial values"<<std::endl;
    sampler.setInitialWalkerPos(initVals, auxVals);
    
    std::cout<<"Running ensemble sampler for autoregressive model"<<std::endl;
    sampler.runMCMC(numSteps);
    
    std::cout<<"Calculating integrated autocorrelation times for autoregressive model"<<std::endl;
    Analysis::AutoCorrCalc<float> acCalc(numParams, numWalkers);
    auto startItt = sampler.getStepIttBegin();
    auto endItt = sampler.getStepIttEnd();
    acCalc.allAutoCorrTime(startItt, endItt, numSteps);
    
    std::cout<<"P0 Calculated AutoCorrelation Time: "<<acCalc.retrieveAutoCorrelationTime(0)<<" Actual:   9"<<std::endl;
    std::cout<<"P1 Calculated AutoCorrelation Time: "<<acCalc.retrieveAutoCorrelationTime(1)<<" Actual:  20"<<std::endl;
    std::cout<<"P2 Calculated AutoCorrelation Time: "<<acCalc.retrieveAutoCorrelationTime(2)<<" Actual:  30"<<std::endl;
    std::cout<<"P3 Calculated AutoCorrelation Time: "<<acCalc.retrieveAutoCorrelationTime(3)<<" Actual:  60"<<std::endl;
    std::cout<<"P4 Calculated AutoCorrelation Time: "<<acCalc.retrieveAutoCorrelationTime(4)<<" Actual: 200"<<std::endl;
}
