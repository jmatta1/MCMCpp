#include<iostream>
#include<sstream>
#include"Utility/NoAction.h"
#include"Movers/Diagnostic/SequenceMove.h"
#include"ParallelEnsembleSampler.h"
using MCMC::ParallelEnsembleSampler;
namespace Mover=MCMC::Mover;

int main(int argc, char* argv[])
{
    if(argc != 3)
    {
        std::cout<<"Usage:\n\t"<<argv[0]<<" <thread_count> <LockType>"<<std::endl;
        std::cout<<"Where LockType can be 0 or 1, 0 for spin locks, 1 for traditional locks"<<std::endl;
        return 1;
    }
    std::istringstream converter;
    converter.str(argv[1]);
    int numThreads = 1;
    converter >> numThreads;
    converter.clear();
    converter.str(argv[2]);
    int option = 0;
    converter >> option;
    const int runNumber = 0;
    const int numWalkers = 2400; //Chosen so that 1, 2, 3, or 4 threads will not have false sharing
    const int numParams = 4;
    const int numSteps = 20000;
    double stepSize[numParams] = {1.0, 2.0, 3.0, 4.0};
    
    std::cout<<"Building initial mover"<<std::endl;
    Mover::SequenceMove<double> mover(numParams, stepSize);
    std::cout<<"Building sampler with: "<<numThreads<<" Threads and "<<numWalkers<<" Walkers"<<std::endl;
    
    bool status = false;
    
    if(option == 0)
    {
        std::cout<<"Sampler will use spin locks"<<std::endl;
        //give the sampler a size larger than normal in order to accomodate the full chain size
        ParallelEnsembleSampler<double, Mover::SequenceMove<double>, MCMC::Utility::NoAction<double>, true>
                sampler(runNumber, numThreads, numWalkers, numParams, mover, 3300000000);
        
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
        status = sampler.runMCMC(numSteps);
        
        auto end = sampler.getStepIttEnd();
        for(auto itt = sampler.getStepIttBegin(); itt != end; ++itt)
        {
            std::cout<<(*itt)[0]<<", "<<(*itt)[1]<<", "<<(*itt)[2]<<", "<<(*itt)[3]<<std::endl;
        }
    }
    else if(option == 1)
    {
        std::cout<<"Sampler will use normal locks"<<std::endl;
        //give the sampler a size larger than normal in order to accomodate the full chain size
        ParallelEnsembleSampler<double, Mover::SequenceMove<double>, MCMC::Utility::NoAction<double>, false>
                sampler(runNumber, numThreads, numWalkers, numParams, mover, 3300000000);
        
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
        status = sampler.runMCMC(numSteps);
        auto end = sampler.getStepIttEnd();
        for(auto itt = sampler.getStepIttBegin(); itt != end; ++itt)
        {
            std::cout<<(*itt)[0]<<", "<<(*itt)[1]<<", "<<(*itt)[2]<<", "<<(*itt)[3]<<std::endl;
        }
    }
    if(status)
    {
        std::cout<<"Sampling completed normally"<<std::endl;
    }
    else
    {
        std::cout<<"Sampling finished when chain ran out of space"<<std::endl;
    }
    return 0;
}
