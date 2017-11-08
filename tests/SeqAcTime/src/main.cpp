#include<iostream>
#include"Movers/Diagnostic/AutoRegressiveMove.h"
#include"Movers/Diagnostic/SequenceMove.h"
#include"Analysis/AutoCorrCalc.h"
#include"EnsembleSampler.h"
#include<fstream>
using std::ofstream;
using MarkovChainMonteCarlo::EnsembleSampler;
namespace Mover=MarkovChainMonteCarlo::Mover;
namespace Analysis=MarkovChainMonteCarlo::Analysis;

static const char* FileNames[5]={"p0.csv", "p1.csv", "p2.csv", "p3.csv", "p4.csv"};

int main()
{
    const int runNumber = 0;
    const int numWalkers = 100;
    const int numParams = 5;
    const int numSteps = 199999;
    double offsets[numParams] = {1.0f, 1.0f, 1.0f, 1.0f, 1.0f};
    double phis[numParams] = {0.8f,  // AC = 9
                             0.904761904762f,  // AC = 20
                             0.9354838709677f,  // AC = 30
                             0.9672131147541f,  // AC = 60
                             0.990050200903734685f};  // AC = 200
    double vars[numParams] = {1.0f, 1.0f, 1.0f, 1.0f, 1.0f};
    
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
    delete[] initVals;
    delete[] auxVals;
}
