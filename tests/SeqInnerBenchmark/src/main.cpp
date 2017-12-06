#include<iostream>
#include"Movers/Diagnostic/SequenceMove.h"
#include"EnsembleSampler.h"
using MCMC::EnsembleSampler;
namespace Mover=MCMC::Mover;

int main()
{
    const int runNumber = 0;
    const int numWalkers = 100;
    const int numParams = 4;
    const int numSteps = 1000000;
    double stepSize[numParams] = {1.0, 2.0, 3.0, 4.0};
    
    std::cout<<"Building initial mover"<<std::endl;
    Mover::SequenceMove<double> mover(numParams, stepSize);
    
    std::cout<<"Setting initial values"<<std::endl;
    double* initVals = new double[numWalkers*numParams];
    double* auxVals = new double[numWalkers];
    for(int i=0; i<numWalkers; ++i)
    {
        auxVals[i] = 0.0;
        for(int j=0; j<numParams; ++j)
        {
            initVals[i*numParams+j] = 0.0;
        }
    }
    
    std::cout<<"Building sampler"<<std::endl;
    //give the sampler a size larger than normal in order to accomodate the full chain size
    EnsembleSampler<double, Mover::SequenceMove<double> > sampler(runNumber, numWalkers, numParams, mover, 3300000000);
    
    std::cout<<"Setting initial values"<<std::endl;
    sampler.setInitialWalkerPos(initVals, auxVals);
    std::cout<<"Deleting initial value arrays"<<std::endl;
    delete[] initVals;
    delete[] auxVals;
    
    std::cout<<"Running ensemble sampler for benchmarking"<<std::endl;
    bool status = sampler.runMCMC(numSteps);
    if(status)
    {
        std::cout<<"Sampling completed normally"<<std::endl;
    }
    else
    {
        std::cout<<"Sampling finished when chain ran out of space"<<std::endl;
    }
}
