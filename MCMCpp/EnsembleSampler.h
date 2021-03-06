/*!*****************************************************************************
********************************************************************************
**
** @copyright Copyright (C) 2017-2018 James Till Matta
**
** This program is distributed in the hope that it will be useful,
** but WITHOUT ANY WARRANTY; without even the implied warranty of
** MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.
** 
********************************************************************************
*******************************************************************************/
#ifndef MCMCPP_ENSEMBLESAMPLER_H
#define MCMCPP_ENSEMBLESAMPLER_H
// includes for C system headers
// includes for C++ system headers
#include<cassert>
// includes from other libraries
// includes from MCMCpp
#include"Chain/Chain.h"
#include"Walker/Walker.h"
#include"Utility/NoAction.h"
#include"Utility/UserOjbectsTest.h"

namespace MCMC
{
/*!
 * @class EnsembleSampler
 * @ingroup Primary
 * @brief The user interface to perform a sequential MCMC sampling
 * @author James Till Matta
 * 
 * @tparam ParamType The floating point type for which this calculation will be carried out
 * @tparam Mover The mover class, current options are StretchMove and WalkMove, StretchMove is
 * computationally much faster, but in some cases it can result in longer auto-correllation times,
 * WalkMove is much more computationally intensive, but it can result in shorter auto-correllation times in some cases
 * @tparam PostStepAction A functor that is called and given a pointer to the chain at the end of every step,
 * it need not be concurrent, it is called by a single thread, though that thread may be different for each call
 */
template<class ParamType, class Mover, class PostStepAction=Utility::NoAction<ParamType> >
class EnsembleSampler
{
public:
    //define some useful typenames
    typedef Chain::Chain<ParamType> ChainType;
    typedef Walker::Walker<ParamType> WalkerType;
    typedef Chain::ChainPsetIterator<ParamType> PsetItt;
    typedef Chain::ChainStepIterator<ParamType> StepItt;
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
     * \param randSeed The random seed used to seed the random number generator (for parallel generators, the stream number will be used)
     * \param numWalker The number of walkers to use, must be a multiple of 2 and exceed 2*numParameter
     * \param numParameter The number of parameters in the problem
     * \param maxChainSizeBytes The maximum size of the sample chain in bytes
     * \param stepAction An instance of the post step action class
     */
    EnsembleSampler(int randSeed, int numWalker, int numParameter, const Mover& move,
                    unsigned long long maxChainSizeBytes=2147483648, PostStepAction* stepAct=nullptr);
    
    /*!
     * @brief ~EnsembleSampler Delete the walker lists then allow the rest to die normally
     */
    ~EnsembleSampler(){delete[] walkerRedSet; delete[] walkerBlkSet;}
    
    /*!
     * \brief Delete copy constructor
     */
    EnsembleSampler(const EnsembleSampler<ParamType, Mover, PostStepAction>& rhs) = delete;
    
    /*!
     * \brief Deleted assignment operator
     */
    EnsembleSampler<ParamType, Mover, PostStepAction>& operator=(const EnsembleSampler<ParamType, Mover, PostStepAction>& rhs) = delete;
    
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
     * \brief runMCMC Runs the Markov chain Monte Carlo for a set number of samples
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
     * \brief getAcceptedSteps Returns the number of steps accepted by the ensemble
     * \return The total number of steps taken by the ensemble
     */
    unsigned long long getAcceptedSteps();
    
    /*!
     * \brief getTotalSteps Returns the total number of steps taken by the ensemble
     * \return The total number of steps taken by the ensembled
     */
    unsigned long long getTotalSteps();
    
    /*!
     * \brief setSlicingMode Is used to change the sampling mode, either using subsampling, or taking all samples
     * \param useSlicing Boolean for whether or not to use subsampling, default is false
     * \param slicingInterval The interval for subsampling
     */
    void setSlicingMode(bool useSlicing=false, int slicingInterval=1);
    
    /*!
     * \brief sliceAndBurnChain Is used to remove burn in samples and select subsamples form the existing chain
     * \param slicingInterval The sub sampling interval to use
     * \param burnIn The number of samples to remove from the beginning of the chain
     */
    void sliceAndBurnChain(int slicingInterval, int burnIn);
    
    /*!
     * \brief getStoredSteps retrieves the number of steps that are stored in the chain
     * \return The number of steps stored in the chain
     */
    int getStoredSteps(){return storedSteps;}
    
    /*!
     * \brief getParamSetIttBegin Gets an iterator pointing to the beginning of the chain
     * that traverses the chain individual parameter set by individual parameter set. Incrementing
     * the iterator gets you the parameter set of the next walker, not paying attention to which parameter set
     * belongs to which step
     * \return A Parameter Set iterator pointed to the beginning of the chain.
     */
    PsetItt getParamSetIttBegin(){return markovChain.getPsetIteratorBegin();}
    /*!
     * \brief getParamSetIttEnd Gets a parameter set iterator pointing to the end of the chain
     * \return A parameter set iterator that points just past the end of the chain
     */
    PsetItt getParamSetIttEnd(){return markovChain.getPsetIteratorEnd();}
    
    /*!
     * \brief getStepIttBegin Gets a step iterator pointing to the chain step by step.
     * Incrementing the iterator takes you to the first walker of the next step
     * \return A Step iterator pointed to the beginning of the chain
     */
    StepItt getStepIttBegin(){return markovChain.getStepIteratorBegin();}
    /*!
     * \brief getStepIttEnd Gets a step iterator that points to the end of chain
     * \return A Step iterator pointed to the end of the chain
     */
    StepItt getStepIttEnd(){return markovChain.getStepIteratorEnd();}
private:
    /*!
     * \brief performStep Performs a step on all of the walkers
     * \param save Whether or not to save the walker's current point, before deciding to update
     */
    void performStep(bool save);
    
    WalkerType* walkerRedSet=nullptr; ///<Set one of the walkers, the sequential mode does not need two sets of walkers, but it is more convenient
    WalkerType* walkerBlkSet=nullptr; ///<Set two of the walkers, the sequential mode does not need two sets of walkers, but it is more convenient
    PostStepAction* stepAction; ///<Action to perform at the end of every step
    ChainType markovChain; ///<The Markov Chain storage class
    Mover moveProposer; ///<The class that proposes a new move position, expects a single numParams constructor parameter
    
    int numParams; ///<The number of parameters being searched on
    int numWalkers; ///<The number of walkers, must be a multiple of 2, and greater than 2*numParams
    int walkersPerSet; ///<The number of walkers in the two sets (numWalkers/2)
    int subSamplingInterval=1; ///<The interval to store samples on if we are subsampling
    int storedSteps=0; ///<The number of steps stored in the chain

    bool subSampling = false; ///<Toggle for performing subsampling
};

template<class ParamType, class Mover, class PostStepAction>
EnsembleSampler<ParamType, Mover, PostStepAction>::
EnsembleSampler(int randSeed, int numWalker, int numParameter, const Mover& move,
                unsigned long long maxChainSizeBytes, PostStepAction* stepAct):
    stepAction(stepAct), markovChain(numWalker, numParameter, maxChainSizeBytes),
    moveProposer(move), numParams(numParameter), numWalkers(numWalker),
    walkersPerSet(numWalker/2)
{
    assert(numWalkers%2 == 0);
    assert(numWalkers > (2*numParams));
    walkerRedSet = new WalkerType[walkersPerSet];
    walkerBlkSet = new WalkerType[walkersPerSet];
    for(int i=0; i<walkersPerSet; ++i)
    {
        walkerRedSet[i].init(&markovChain, i,   numParams);
        walkerBlkSet[i].init(&markovChain, i+walkersPerSet, numParams);
    }
    //set the random seed to what was given to the user and the stream number to 0
    moveProposer.setPrng(randSeed, 0);
}

template<class ParamType, class Mover, class PostStepAction>
void EnsembleSampler<ParamType, Mover, PostStepAction>::setInitialWalkerPos(ParamType* positions, ParamType* auxValues)
{
    for(int i=0; i<walkersPerSet; ++i)
    {
        walkerRedSet[i].setFirstPoint(positions+(i*numParams), auxValues[i], true);
        walkerBlkSet[i].setFirstPoint(positions+((i+walkersPerSet)*numParams), auxValues[(i+walkersPerSet)], true);
    }
    ++storedSteps;
    markovChain.incrementChainStep();
}

template<class ParamType, class Mover, class PostStepAction>
void EnsembleSampler<ParamType, Mover, PostStepAction>::storeCurrentWalkerPositions()
{
    for(int i=0; i<walkersPerSet; ++i)
    {
        walkerRedSet[i].storeCurrentPoint();
        walkerBlkSet[i].storeCurrentPoint();
    }
    ++storedSteps;
    markovChain.incrementChainStep();
}

template<class ParamType, class Mover, class PostStepAction>
ParamType EnsembleSampler<ParamType, Mover, PostStepAction>::getAcceptanceFraction()
{
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
unsigned long long EnsembleSampler<ParamType, Mover, PostStepAction>::getAcceptedSteps()
{
    unsigned long long accepted = 0ULL;
    for(int i=0; i<walkersPerSet; ++i)
    {
        accepted += walkerRedSet[i].getAcceptedProposals();
        accepted += walkerBlkSet[i].getAcceptedProposals();
    }
    return accepted;
}

template<class ParamType, class Mover, class PostStepAction>
unsigned long long EnsembleSampler<ParamType, Mover, PostStepAction>::getTotalSteps()
{
    unsigned long long total = 0ULL;
    for(int i=0; i<walkersPerSet; ++i)
    {
        total += walkerRedSet[i].getTotalSteps();
        total += walkerBlkSet[i].getTotalSteps();
    }
    return total;
}

template<class ParamType, class Mover, class PostStepAction>
bool EnsembleSampler<ParamType, Mover, PostStepAction>::runMCMC(int numSteps)
{
    if(!subSampling)
    {
        for(int i=0; i<numSteps; ++i)
        {
            performStep(true);
            ++storedSteps;
            if(Chain::IncrementStatus::EndOfChain == markovChain.incrementChainStep()) return false;
        }
    }
    else
    {
        for(int i=0; i<numSteps; ++i)
        {//since the initialization is counted do the subSamplingInterval-1 steps without storing
            for(int j=1; j<subSamplingInterval; ++j)
            {
                performStep(false);
            }
            performStep(true);
            ++storedSteps;
            if(Chain::IncrementStatus::EndOfChain == markovChain.incrementChainStep()) return false;
        }
    }
    return true;
}

template<class ParamType, class Mover, class PostStepAction>
void EnsembleSampler<ParamType, Mover, PostStepAction>::reset()
{
    storedSteps = 0;
    markovChain.resetChain();
    for(int i=0; i<walkersPerSet; ++i)
    {
        walkerRedSet[i].resetSteps();
        walkerBlkSet[i].resetSteps();
    }
}

template<class ParamType, class Mover, class PostStepAction>
void EnsembleSampler<ParamType, Mover, PostStepAction>::setSlicingMode(bool useSlicing, int slicingInterval)
{
    assert(slicingInterval>0);
    subSampling=useSlicing;
    subSamplingInterval = slicingInterval;
}

template<class ParamType, class Mover, class PostStepAction>
void EnsembleSampler<ParamType, Mover, PostStepAction>::sliceAndBurnChain(int slicingInterval, int burnIn)
{
    assert(slicingInterval>0);
    assert(burnIn>=0);
    markovChain.resetChainForSubSampling(burnIn, slicingInterval);
    storedSteps = markovChain.getStoredStepCount();
}

template<class ParamType, class Mover, class PostStepAction>
void EnsembleSampler<ParamType, Mover, PostStepAction>::performStep(bool save)
{
    //first update all the red set from the black set
    for(int i=0; i<walkersPerSet; ++i)
    {
        moveProposer.updateWalker(walkerRedSet[i], walkerBlkSet, walkersPerSet, save);
    }
    //now update all the black set from the red set
    for(int i=0; i<walkersPerSet; ++i)
    {
        //first get a proposed new point to move to
        moveProposer.updateWalker(walkerBlkSet[i], walkerRedSet, walkersPerSet, save);
    }
    //now perform the poststep action
    if(stepAction != nullptr)
    {
        stepAction->performAction(markovChain.getStepIteratorBegin(), markovChain.getStepIteratorEnd());
    }
}

}
#endif  //MCMCPP_ENSEMBLESAMPLER_H
