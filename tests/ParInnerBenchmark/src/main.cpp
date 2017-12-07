#include<iostream>
#include<sstream>
#include"Movers/Diagnostic/SequenceMove.h"
#include"ParallelEnsembleSampler.h"
using MCMC::ParallelEnsembleSampler;
namespace Mover=MCMC::Mover;

int main(int argc, char* argv[])
{
    if(argc != 2)
    {
        std::cout<<"Usage:\n\t"<<argv[0]<<" <thread_count>"<<std::endl;
        return 1;
    }
    std::istringstream converter;
    converter.str(argv[1]);
    int numThreads = 1;
    converter >> numThreads;
    const int runNumber = 0;
    const int numWalkers = 240; //Chosen so that 1, 2, 3, or 4 threads will not have false sharing
    const int numParams = 4;
    const int numSteps = 200000;
    double stepSize[numParams] = {1.0, 2.0, 3.0, 4.0};
    
    std::cout<<"Building initial mover"<<std::endl;
    Mover::SequenceMove<double> mover(numParams, stepSize);
    
    std::cout<<"Building sampler with: "<<numThreads<<" Threads and "<<numWalkers<<" Walkers"<<std::endl;
    //give the sampler a size larger than normal in order to accomodate the full chain size
    ParallelEnsembleSampler<double, Mover::SequenceMove<double> > sampler(runNumber, numThreads, numWalkers, numParams, mover, 3300000000);
    
    std::cout<<"Getting initial values"<<std::endl;
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
    
    std::cout<<"Setting initial values"<<std::endl;
    sampler.setInitialWalkerPos(initVals, auxVals);
    std::cout<<"Deleting initial value arrays"<<std::endl;
    delete[] initVals;
    delete[] auxVals;
    
    std::cout<<"Running parallel ensemble sampler for benchmarking"<<std::endl;
    bool status = sampler.runMCMC(numSteps);
    
    if(status)
    {
        std::cout<<"Sampling completed normally"<<std::endl;
    }
    else
    {
        std::cout<<"Sampling finished when chain ran out of space"<<std::endl;
    }
    
    /*auto end = sampler.getStepIttEnd();
    for(auto itt = sampler.getStepIttBegin(); itt != end; ++itt)
    {
        std::cout<<(*itt)[0]<<", "<<(*itt)[1]<<", "<<(*itt)[2]<<", "<<(*itt)[3]<<std::endl;
    }*/
}
