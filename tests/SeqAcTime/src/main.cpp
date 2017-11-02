#include<iostream>
#include"Movers/AutoRegressiveMover.h"
#include"Analysis/AutoCorrCalc.h"
#include"EnsembleSampler.h"
using MarkovChainMonteCarlo::EnsembleSampler;
namespace Mover=MarkovChainMonteCarlo::Mover;
namespace Analysis=MarkovChainMonteCarlo::Analysis;

int main()
{
    float offsets[5] = {1.0f, 1.0f, 1.0f, 1.0f, 1.0f};
    float phis[5] = {0.8f,  // AC = 9
                     0.904761904762f,  // AC = 20
                     0.9354838709677f,  // AC = 30
                     0.9672131147541f,  // AC = 60
                     0.990050200903734685f};  // AC = 200
    float vars[5] = {1.0f, 1.0f, 1.0f, 1.0f, 1.0f};
    Mover::AutoRegressiveMove<float> mover(5, offsets, phis, vars);
    EnsembleSampler<float, Mover::AutoRegressiveMove<float> > sampler(0, 100, 5, mover);
    std::cout<<"Running ensemble sampler for autoregressive model"<<std::endl;
    sampler.runMCMC(200000);
    auto startItt = sampler.getStepIttBegin();
    auto endItt = sampler.getStepIttEnd();
    std::cout<<"Calculating integrated autocorrelation times for autoregressive model"<<std::endl;
    Analysis::AutoCorrCalc<float> acCalc(5, 100);
    acCalc.allAutoCorrTime(startItt, endItt, 200000);
    std::cout<<"Paramater 0 Calculated AutoCorrelation Time: "<<acCalc.retrieveAutoCorrelationTime(0)<<" Actual: 9"<<std::endl;
    std::cout<<"Paramater 1 Calculated AutoCorrelation Time: "<<acCalc.retrieveAutoCorrelationTime(1)<<" vs. 20"<<std::endl;
    std::cout<<"Paramater 2 Calculated AutoCorrelation Time: "<<acCalc.retrieveAutoCorrelationTime(2)<<" vs. 30"<<std::endl;
    std::cout<<"Paramater 3 Calculated AutoCorrelation Time: "<<acCalc.retrieveAutoCorrelationTime(3)<<" vs. 60"<<std::endl;
    std::cout<<"Paramater 4 Calculated AutoCorrelation Time: "<<acCalc.retrieveAutoCorrelationTime(4)<<" vs. 200"<<std::endl;
}
