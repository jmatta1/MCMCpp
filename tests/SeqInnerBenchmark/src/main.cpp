#include<iostream>
#include"Movers/Diagnostic/SequenceMove.h"
#include"EnsembleSampler.h"
using MarkovChainMonteCarlo::EnsembleSampler;
namespace Mover=MarkovChainMonteCarlo::Mover;

int main()
{
    const int runNumber = 0;
    const int numWalkers = 100;
    const int numParams = 5;
    const int numSteps = 1000000;
    double stepSize[numParams] = {1.0, 2.0, 3.0, 4.0, 5.0};
    
    std::cout<<"Building initial mover"<<std::endl;
    Mover::SequenceMove<double> mover(numParams, stepSize);
    
    std::cout<<"Setting initial values"<<std::endl;
    double* initVals = new double[numWalkers*numParams];
    double* auxVals = new double[numWalkers];
    for(int i=0; i<numWalkers; ++i)
    {
        auxVals[i] = 0.0;
        int limit = ((i*numParams) + numParams);
        for(int j=(i*numParams); j<limit; ++j)
        {
            initVals[j] = 0.0;
        }
    }
    
    std::cout<<"Building sampler"<<std::endl;
    EnsembleSampler<double, Mover::SequenceMove<double> > sampler(runNumber, numWalkers, numParams, mover);
    
    std::cout<<"Setting initial values"<<std::endl;
    sampler.setInitialWalkerPos(initVals, auxVals);
    std::cout<<"Deleting initial value arrays"<<std::endl;
    delete[] initVals;
    delete[] auxVals;
    
    std::cout<<"Running ensemble sampler for benchmarking"<<std::endl;
    sampler.runMCMC(numSteps);
}
