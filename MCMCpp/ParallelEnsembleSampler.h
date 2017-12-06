/*!*****************************************************************************
********************************************************************************
**
** @copyright Copyright (C) 2017 James Till Matta
**
** This program is distributed in the hope that it will be useful,
** but WITHOUT ANY WARRANTY; without even the implied warranty of
** MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.
** 
********************************************************************************
*******************************************************************************/
#ifndef MCMC_PARALLELENSEMBLESAMPLER_H
#define MCMC_PARALLELENSEMBLESAMPLER_H
// includes for C system headers
// includes for C++ system headers
#include<cassert>
#include<mutex>
#include<thread>
// includes from other libraries
// includes from MCMC
#include"Chain/Chain.h"
#include"Walker/Walker.h"
#include"Utility/NoAction.h"
#include"Utility/UserOjbectsTest.h"
#include"Threading/RedBlkCtrler.h"
#include"Threading/RedBlkUpdater.h"
#include"Threading/ThreadWrapper.h"

namespace MCMC
{

/*!
 * @class ParallelEnsembleSampler
 * @ingroup Primary
 * @brief The user interface to perform a parallel MCMC sampling
 * @author James Till Matta
 * 
 * @tparam ParamType The floating point type for which this calculation will be
 * carried out
 * @tparam Mover The mover class, current options are StretchMove and WalkMove,
 * StretchMove is computationally much faster, but in some cases it can result
 * in longer auto-correllation times, WalkMove is much more computationally
 * intensive, but it can result in shorter auto-correllation times in some cases
 * The mover class will need to be reentrant, that is to say that any shared
 * state between copies *must* be read only. This allows multiple threads to
 * use multiple copies while still maintaining low communication
 * @tparam PostStepAction A functor that is called and given an iterattor to the
 * chain start and end at the end of every step, it need not be concurrent or
 * reentrant, it is only ever called by a single thread, though the thread 
 * making the call may be different between calls
 * 
 * **Recommendation**: numWalkers should be chosen such that
 * (numWalkers/(2*numThreads))*numParam*sizeof(ParamType) is an integer
 * multiple of 128 bytes. This ensures that individual threads are not sharing
 * data that is split across cachelines, preventing false sharing of data
 * between threads, allowing a faster calculation (because the cache is not
 * being invalidated for one thread when another thread updates the chain it
 * is working on that happens to share a cache line with the chain of the
 * first thread, invalidating the cache even though neither thread is
 * examining the same data)
 */
template<class ParamType, class Mover, class PostStepAction=Utility::NoAction<ParamType> >
class ParallelEnsembleSampler
{
public:
    //define some useful typenames
    typedef Chain::Chain<ParamType> ChainType;
    typedef Walker::Walker<ParamType> WalkerType;
    typedef Threading::RedBlackCtrler<ParamType, PostStepAction> CtrlType;
    typedef Chain::ChainPsetIterator<ParamType> PsetItt;
    typedef Chain::ChainStepIterator<ParamType> StepItt;
    typedef Threading::RedBlkUpdater<ParamType, Mover, PostStepAction> ThreadObjectType;

    //perform static checks of the users classes to ensure that they have the needed member functions for their role
    static_assert(Utility::CheckCalcUpdateWalker<Mover, void, WalkerType&, WalkerType*, int, bool>(),
                  "The Mover class does not have the necessary member function with signature:\n"
                  "  'void updateWalker(Walker::Walker<ParamType>&, Walker::Walker<ParamType>*, int, bool)'");
    static_assert(Utility::CheckPerformAction<PostStepAction, void, const StepItt&, const StepItt&>(),
                  "The PostStepAction class does not have the necessary member function with signature:\n"
                  "  'void PerformAction(Chain::ChainPsetIterator<ParamType>& start, Chain::ChainPsetIterator<ParamType>& end)'");
    static_assert(std::is_copy_constructible<Mover>::value, "The Mover class needs to be copy constructible.");
    static_assert(std::is_copy_constructible<PostStepAction>::value, "The PostStepAction class needs to be copy constructible.");
    
    /*!
     * \brief EnsembleSampler Constructs the ensemble sampler
     * \param randSeed The random seed used to seed the random number generator (for parallel generators, the stream number will be used to differentiate)
     * \param numThreads The number of threads to utilize to speed the calculation
     * \param numWalker The number of walkers to use, must be a multiple of 2 and exceed 2*numParameter.
     * \param numParameter The number of parameters in the problem
     * \param maxChainSizeBytes The maximum size of the sample chain in bytes
     * \param stepAction An instance of the post step action class
     * 
     * **Recommendation**: numWalkers should be chosen such that
     * (numWalkers/(2*numThreads))*numParam*sizeof(ParamType) is an integer
     * multiple of 128 bytes. This ensures that individual threads are not sharing
     * data that is split across cachelines, preventing false sharing of data
     * between threads, allowing a faster calculation (because the cache is not
     * being invalidated for one thread when another thread updates the chain it
     * is working on that happens to share a cache line with the chain of the
     * first thread, invalidating the cache even though neither thread is
     * examining the same data)
     */
    ParallelEnsembleSampler(int randSeed, int threadCount, int numWalker, int numParameter, const Mover& move,
                            unsigned long long maxChainSizeBytes=2147483648, PostStepAction* stepAct=nullptr);
    
    /*!
     * @brief ~ParallelEnsembleSampler Delete the walker lists and temp parameter set then allow the rest to die normally
     */
    ~ParallelEnsembleSampler();
    
    /*!
     * \brief Delete copy constructor
     */
    ParallelEnsembleSampler(const ParallelEnsembleSampler<ParamType, Mover, PostStepAction>& rhs) = delete;
    
    /*!
     * \brief Deleted assignment operator
     */
    ParallelEnsembleSampler<ParamType, Mover, PostStepAction>& operator=(const ParallelEnsembleSampler<ParamType, Mover, PostStepAction>& rhs) = delete;
    
    /*!
     * \brief setInitialWalkerPos Gives an initial position to every walker.
     * \param positions An array of floating point types with length numParameters*numWalker representing a starting point for every walker
     * \param auxValues An array of floating point types with length numWalker representing the auxillary data for every walker
     */
    void setInitialWalkerPos(ParamType* positions, ParamType* auxValues);
    
    /*!
     * \brief storeWalkerPositions Takes the current position for each walker and stores it into the chain
     * 
     * This is useful for if you have reset the chain but wish to use the current
     * position that each walker is at
     */
    void storeCurrentWalkerPositions();
    
    /*!
     * \brief runMCMC Runs the Markov chain Monte Carlo for a set number of samples across multiple threads
     * \param numSteps The number of steps to store
     * \return True if the sampling went to the end before the max chain size was reached, false otherwise
     * 
     * In normal sampling mode, this function will run the ensemble for numSamples per walker
     * In subSampling mode, this function will run the ensemble numSamples*subSamplingInterval and store numSamples
     */
    bool runMCMC(int numSteps);
    
    /*!
     * \brief reset Returns the sampler to it's original state except that the walkers retain their current positions
     */
    void reset();
    
    /*!
     * \brief getAcceptanceFraction Calculated the acceptance fraction for the ensemble over the current steps
     * \return The fraction of steps where the proposal was accepted
     */
    ParamType getAcceptanceFraction();
    
    /*!
     * \brief setSamplingMode Is used to change the sampling mode, either using subsampling, or taking all samples
     * \param subSamplingInt The interval for subsampling, if greater than 1, then subSamplingInt-1 samples will be skipped for every 1 sample taken
     * \param burnIn The number of points at the beginning to discard for burnin
     */
    void setSamplingMode(int subSamplingInt=1, int burnIn=0);
    
    /*!
     * \brief getStoredSteps retrieves the number of steps that are stored in the chain
     * \return The number of steps stored in the chain
     */
    int getStoredSteps(){std::unique_lock<std::mutex> lock(samplerMutex); return markovChain.getStoredStepCount();}
    
    /*!
     * \brief getParamSetIttBegin Gets an iterator pointing to the beginning of the chain
     * that traverses the chain individual parameter set by individual parameter set. Incrementing
     * the iterator gets you the parameter set of the next walker, not paying attention to which parameter set
     * belongs to which step
     * \return A Parameter Set iterator pointed to the beginning of the chain.
     */
    PsetItt getParamSetIttBegin(){std::unique_lock<std::mutex> lock(samplerMutex); return markovChain.getPsetIteratorBegin();}
    /*!
     * \brief getParamSetIttEnd Gets a parameter set iterator pointing to the end of the chain
     * \return A parameter set iterator that points just past the end of the chain
     */
    PsetItt getParamSetIttEnd(){std::unique_lock<std::mutex> lock(samplerMutex); return markovChain.getPsetIteratorEnd();}
    
    /*!
     * \brief getStepIttBegin Gets a step iterator pointing to the chain step by step.
     * Incrementing the iterator takes you to the first walker of the next step
     * \return A Step iterator pointed to the beginning of the chain
     */
    StepItt getStepIttBegin(){std::unique_lock<std::mutex> lock(samplerMutex); return markovChain.getStepIteratorBegin();}
    /*!
     * \brief getStepIttEnd Gets a step iterator that points to the end of chain
     * \return A Step iterator pointed to the end of the chain
     */
    StepItt getStepIttEnd(){std::unique_lock<std::mutex> lock(samplerMutex); return markovChain.getStepIteratorEnd();}
private:
    ChainType markovChain; ///<The Markov Chain storage class
    CtrlType controller; ///<Object used to control threads
    std::mutex samplerMutex; ///<Used to prevent calls during sampling, i.e. only one thread can access this object at a time
    WalkerType* walkerRedSet=nullptr; ///<Set one of the walkers
    WalkerType* walkerBlkSet=nullptr; ///<Set two of the walkers
    ThreadObjectType** threadObjects=nullptr; ///<Set of callable objects used to construct the threads
    std::thread** threadGroup=nullptr; ///<Set of thread objects which reprresent the threads of execution
    
    
    int numThreads; ///<The number of threads to spread the work across
    int numParams; ///<The number of parameters being searched on
    int numWalkers; ///<The number of walkers, must be a multiple of 2, and greater than 2*numParams
    int walkersPerSet; ///<The number of walkers in each of the two sets (numWalkers/2)
    int subSamplingInterval=1; ///<The interval to store samples on if we are subsampling
};

template<class ParamType, class Mover, class PostStepAction>
ParallelEnsembleSampler<ParamType, Mover, PostStepAction>::ParallelEnsembleSampler(
        int randSeed, int threadCount, int numWalker, int numParameter, const Mover& move,
        unsigned long long maxChainSizeBytes, PostStepAction* stepAct):
    markovChain(numWalker, numParameter, maxChainSizeBytes),
    controller(threadCount, markovChain, stepAct), samplerMutex(), numThreads(threadCount),
    numParams(numParameter), numWalkers(numWalker), walkersPerSet(numWalker/2)
{
    std::unique_lock<std::mutex> lock(samplerMutex);
    assert(numWalkers%2 == 0);
    assert(numWalkers > (2*numParams));
    //set up the walkers
    walkerRedSet = new WalkerType[walkersPerSet];
    walkerBlkSet = new WalkerType[walkersPerSet];
    for(int i=0; i<walkersPerSet; ++i)
    {
        walkerRedSet[i].init(&markovChain, i,   numParams);
        walkerBlkSet[i].init(&markovChain, i+walkersPerSet, numParams);
    }
    //now set up the threads
    threadObjects = new ThreadObjectType*[numThreads];
    threadGroup = new std::thread*[numThreads];
    int offset = 0;
    int nominalWalkersPerThread = (walkersPerSet/numThreads);
    int excess = (walkersPerSet%numThreads);
    int size = ((excess == 0) ? nominalWalkersPerThread : nominalWalkersPerThread+1);
    for(int i=0; i<numThreads; ++i)
    {
        //check if we have finished removing the excess
        if(i==excess) size = nominalWalkersPerThread;
        typename ThreadObjectType::WalkerInfo walkerSets = std::make_tuple(walkerRedSet, walkerBlkSet, walkersPerSet, walkersPerSet);
        typename ThreadObjectType::WalkerInfo updateSets = std::make_tuple(walkerRedSet+offset, walkerBlkSet+offset, size, size);
        threadObjects[i] = new ThreadObjectType(randSeed, i, walkerSets, updateSets, move, controller);
        //now update the offset
        offset += size;
        //now create the thread
        threadGroup[i] = new std::thread(Threading::ThreadWrapper<ThreadObjectType>(threadObjects[i]));
    }
    //now wait for all the threads to come to the wait state
    while(controller.getNumWorkersWaiting() < numThreads)
    {
        std::this_thread::sleep_for(std::chrono::microseconds(10));
    }
}

template<class ParamType, class Mover, class PostStepAction>
ParallelEnsembleSampler<ParamType, Mover, PostStepAction>::~ParallelEnsembleSampler()
{
    std::unique_lock<std::mutex> lock(samplerMutex);
    //terminate the workers and wait until they all acknowledge the terminate
    controller.terminateWorkers();
    //now that the threads have all acknowledged the terminate, call join on each
    //and then when that returns delete the thread object and the function object
    for(int i=0; i<numThreads; ++i)
    {
        threadGroup[i]->join();
        delete threadGroup[i];
        delete threadObjects[i];
    }
    //delete the arrays needed for various things
    delete[] threadGroup;
    delete[] threadObjects;
    delete[] walkerRedSet;
    delete[] walkerBlkSet;
}

template<class ParamType, class Mover, class PostStepAction>
bool ParallelEnsembleSampler<ParamType, Mover, PostStepAction>::runMCMC(int numSteps)
{
    std::unique_lock<std::mutex> lock(samplerMutex);
    controller.runSampling(numSteps, subSamplingInterval-1);
    return controller.samplingComplete();
}

template<class ParamType, class Mover, class PostStepAction>
ParamType ParallelEnsembleSampler<ParamType, Mover, PostStepAction>::getAcceptanceFraction()
{
    std::unique_lock<std::mutex> lock(samplerMutex);
    unsigned long long accepted = 0ULL;
    unsigned long long total = 0ULL;
    for(int i=0; i<walkersPerSet; ++i)
    {
        accepted += walkerRedSet[i].getAcceptedProposals();
        accepted += walkerBlkSet[i].getAcceptedProposals();
        total += walkerRedSet[i].getTotalSteps();
        total += walkerBlkSet[i].getTotalSteps();
    }
    ParamType fraction = static_cast<ParamType>(accepted)/static_cast<ParamType>(total);
    return fraction;
}

template<class ParamType, class Mover, class PostStepAction>
void ParallelEnsembleSampler<ParamType, Mover, PostStepAction>::reset()
{
    std::unique_lock<std::mutex> lock(samplerMutex); 
    markovChain.resetChain();
    for(int i=0; i<walkersPerSet; ++i)
    {
        walkerRedSet[i].resetSteps();
        walkerBlkSet[i].resetSteps();
    }
}

template<class ParamType, class Mover, class PostStepAction>
void ParallelEnsembleSampler<ParamType, Mover, PostStepAction>::setSamplingMode(int subSamplingInt, int burnIn)
{
    assert(subSamplingInt>0);
    assert(burnIn>=0);
    std::unique_lock<std::mutex> lock(samplerMutex);
    subSamplingInterval = subSamplingInt;
    markovChain.resetChainForSubSampling(burnIn, subSamplingInt);
}

template<class ParamType, class Mover, class PostStepAction>
void ParallelEnsembleSampler<ParamType, Mover, PostStepAction>::setInitialWalkerPos(ParamType* positions, ParamType* auxValues)
{
    std::unique_lock<std::mutex> lock(samplerMutex);
    for(int i=0; i<walkersPerSet; ++i)
    {
        walkerRedSet[i].setFirstPoint(positions+(i*numParams), auxValues[i], true);
        walkerBlkSet[i].setFirstPoint(positions+((i+walkersPerSet)*numParams), auxValues[(i+walkersPerSet)], true);
    }
    markovChain.incrementChainStep();
}

template<class ParamType, class Mover, class PostStepAction>
void ParallelEnsembleSampler<ParamType, Mover, PostStepAction>::storeCurrentWalkerPositions()
{
    std::unique_lock<std::mutex> lock(samplerMutex);
    for(int i=0; i<walkersPerSet; ++i)
    {
        walkerRedSet[i].storeCurrentPoint();
        walkerBlkSet[i].storeCurrentPoint();
    }
    markovChain.incrementChainStep();
}

}
#endif  //MCMC_PARALLELENSEMBLESAMPLER_H
