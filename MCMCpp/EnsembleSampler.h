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
#ifndef MCMC_ENSEMBLESAMPLER_H
#define MCMC_ENSEMBLESAMPLER_H
// includes for C system headers
// includes for C++ system headers
#include<cassert>
// includes from other libraries
// includes from MCMC
#include"Chain/Chain.h"
#include"Walker/Movers/StretchMove.h"
#include"Walker/Movers/WalkMove.h"
#include"Walker/Walker.h"
#include"Utility/NoAction.h"
#include"Utility/GwDistribution.h"
#include"Utility/UserOjbectsTest.h"

namespace MarkovChainMonteCarlo
{
/*!
 * @class EnsembleSampler
 * @ingroup Primary
 * @brief The user interface to perform a sequential MCMC sampling
 * @author James Till Matta
 * 
 * @tparam ParamType The floating point type for which this calculation will be carried out
 * @tparam PostProbCalculator The calculator of the log of the posterior probability function.
 * This can be calculated as the log of the likelihood function summed with the log of the prior-probability.
 * @tparam BlockSize Number of steps to store in a block of the linked list that stores the chain  
 * @tparam PostStepAction A functor that is called and given a pointer to the chain at the end of every step,
 * it need not be concurrent, it is called by a single thread, though that thread may be different for each call
 * @tparam Mover The mover class, current options are StretchMove and WalkMove, StretchMove is
 * computationally much faster, but in some cases it can result in longer auto-correllation times,
 * WalkMove is much more computationally intensive, but it can result in shorter auto-correllation times in some cases
 * @tparam CustomDistribution The class that will be used to set the get sample from custom
 * distribution function in the MultiSampler, currently only affects StretchMove 
 */
template<class ParamType,
         class PostProbCalculator,
         int BlockSize=1000,
         class PostStepAction=Utility::NoAction<ParamType, BlockSize>,
         class CustomDistribution=Utility::GwDistribution<ParamType, 2, 1>,
         class Mover=Walker::StretchMove<ParamType, BlockSize, CustomDistribution, PostProbCalculator> >
class EnsembleSampler
{
public:
    //define some useful typenames
    typedef Chain::Chain<ParamType, BlockSize> ChainType;
    typedef Walker::Walker<ParamType, BlockSize, Mover, PostProbCalculator> WalkerType;
    typedef Utility::MultiSampler<ParamType, CustomDistribution> PrngType;
    typedef Chain::ChainPsetIterator<ParamType, BlockSize> PsetItt;
    typedef Chain::ChainStepIterator<ParamType, BlockSize> StepItt;
    //perform static checks of the users classes to ensure that they have the needed member functions for their role
    static_assert(Utility::CheckPerformAction<PostStepAction, void, PsetItt, PsetItt>(),
                  "The PostStepAction class does not have the necessary member function with signature:\n"
                  "  'void PerformAction(Chain::ChainPsetIterator<ParamType, BlockSize>& start, Chain::ChainPsetIterator<ParamType, BlockSize>& end)'");
    static_assert(Utility::CheckFunctor<CustomDistribution, ParamType, ParamType>(),
                  "The CustomDistribution class does not have the necessary member function with signature:\n"
                  "  'ParamType operator()(ParamType)'");
    static_assert(std::is_trivially_constructible<CustomDistribution>::value, "The CustomDistribution class needs to be trivially constructible.");
    static_assert(std::is_copy_constructible<Mover>::value, "The Mover class needs to be copy constructible.");
    
    /*!
     * \brief EnsembleSampler Constructs the ensemble sampler
     * \param runNumber RunNumber, used to seed the random number generator
     * \param numWalker The number of walkers to use, must be a multiple of 2 and exceed 2*numParameter
     * \param numParameter The number of parameters in the problem
     * \param maxChainSizeBytes The maximum size of the sample chain in bytes
     * \param stepAction An instance of the post step action class
     */
    EnsembleSampler(int runNumber, int numWalker, int numParameter, unsigned long long maxChainSizeBytes, PostStepAction& stepAct);
    
    /*!
     * @brief ~EnsembleSampler Delete the walker lists and temp parameter set then allow the rest to die normally
     */
    ~EnsembleSampler(){delete[] walkerRedSet; delete[] walkerBlkSet; delete[] proposalPoint;}
    
    /*!
     * \brief setInitialWalkerPos Gives an initial position to every walker.
     * \param positions An array of floating point types with length numParameters*numWalker
     */
    void setInitialWalkerPos(ParamType* positions);
    
    /*!
     * \brief runMCMC Runs the Markov chain Monte Carlo for a set number of samples
     * \param numSteps The number of steps to store
     * 
     * In normal sampling mode, this function will run the ensemble for numSamples per walker
     * In subSampling mode, this function will run the ensemble numSamples*subSamplingInterval and store numSamples
     */
    void runMCMC(int numSteps);
    
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
     * \param useSubSampling Boolean for whether or not to use subsampling, default is false
     * \param subSamplingInt The interval for subsampling
     * \param burnIn The number of points at the beginning to discard for burnin
     */
    void setSamplingMode(bool useSubSampling=false, int subSamplingInt=1, int burnIn=0);
    
    /*!
     * \brief calculateAutocorrelationTimes Uses an AutoCorrCalc object to calculate the autocorrelation times for the present chain
     * \return True if the chain was long enough to extract an autocorrelation time with these parameters, False otherwise
     */
    bool calculateAutocorrelationTimes();
    
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
    
    int numParams; ///<The number of parameters being searched on
    int numWalkers; ///<The number of walkers, must be a multiple of 2, and greater than 2*numParams
    int walkersPerSet; ///<The number of walkers in the two sets (numWalkers/2)
    int subSamplingInterval=1; ///<The interval to store samples on if we are subsampling
    int storedSteps=0; ///<The number of steps stored in the chain
    
    ChainType markovChain; ///<The Markov Chain storage class
    WalkerType* walkerRedSet; ///<Set one of the walkers, the sequential mode does not need two sets of walkers, but it is more convenient
    WalkerType* walkerBlkSet; ///<Set two of the walkers, the sequential mode does not need two sets of walkers, but it is more convenient
    ParamType* proposalPoint; ///<Internal parameter set needed to move proposal from mover to walker
    Mover moveProposer; ///<The class that proposes a new move position, expects a single numParams constructor parameter
    PostProbCalculator calc; ///<The class that calculates the posterior probability for a given set of points
    PrngType prng; ///<The pseudorandom number generator used to generate needed random numbers

    bool subSampling = false; ///<Toggle for performing subsampling

    PostStepAction& stepAction; ///<Action to perform at the end of every step
};

template<class ParamType, class PostProbCalculator, int BlockSize,
         class PostStepAction,class CustomDistribution, class Mover>
EnsembleSampler<ParamType, PostProbCalculator, BlockSize, PostStepAction, CustomDistribution, Mover>::
EnsembleSampler(int runNumber, int numWalker, int numParameter, unsigned long long maxChainSizeBytes, PostStepAction& stepAct):
    numParams(numParameter), numWalkers(numWalker), walkersPerSet(numWalker/2),
    markovChain(numWalkers, numParams, maxChainSizeBytes),
    moveProposer(numParams), prng(runNumber), stepAction(stepAct)
{
    assert(numWalkers%2 == 0);
    assert(numWalkers > (2*numParams));
    //allocate the walker arrays
    //for the red set, allocate half the walkers
    walkerRedSet = new WalkerType[walkersPerSet];
    walkerBlkSet = new WalkerType[walkersPerSet];
    proposalPoint = new ParamType[numParams];
    for(int i=0; i<walkersPerSet; ++i)
    {
        walkerRedSet[i].init(&markovChain, 2*i,   numParams);
        walkerBlkSet[i].init(&markovChain, 2*i+1, numParams);
    }
}

template<class ParamType, class PostProbCalculator, int BlockSize,
         class PostStepAction,class CustomDistribution, class Mover>
void EnsembleSampler<ParamType, PostProbCalculator, BlockSize, PostStepAction, CustomDistribution, Mover>::setInitialWalkerPos(ParamType* positions)
{
    for(int i=0; i<walkersPerSet; ++i)
    {
        walkerRedSet[i].setFirstPoint(positions+2*i*numParams, calc);
        walkerBlkSet[i].setFirstPoint(positions+(2*i+1)*numParams, calc);
    }
    ++storedSteps;
    markovChain.incrementChainStep();
}

template<class ParamType, class PostProbCalculator, int BlockSize,
         class PostStepAction,class CustomDistribution, class Mover>
ParamType EnsembleSampler<ParamType, PostProbCalculator, BlockSize, PostStepAction, CustomDistribution, Mover>::getAcceptanceFraction()
{
    unsigned long long accepted;
    unsigned long long total;
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

template<class ParamType, class PostProbCalculator, int BlockSize,
         class PostStepAction,class CustomDistribution, class Mover>
void EnsembleSampler<ParamType, PostProbCalculator, BlockSize, PostStepAction, CustomDistribution, Mover>::runMCMC(int numSteps)
{
    if(!subSampling)
    {
        for(int i=0; i<numSteps; ++i)
        {
            performStep(true);
            ++storedSteps;
            markovChain.incrementChainStep();
        }
    }
    else
    {
        for(int i=0; i<numSteps; ++i)
        {
            performStep(true);
            ++storedSteps;
            markovChain.incrementChainStep();
            for(int j=1; j<subSamplingInterval; ++j)
            {
                performStep(false);
            }
        }
    }
}

template<class ParamType, class PostProbCalculator, int BlockSize,
         class PostStepAction,class CustomDistribution, class Mover>
void EnsembleSampler<ParamType, PostProbCalculator, BlockSize, PostStepAction, CustomDistribution, Mover>::reset()
{
    storedSteps = 0;
    markovChain.resetChain();
    for(int i=0; i<walkersPerSet; ++i)
    {
        walkerRedSet[i].resetSteps();
        walkerBlkSet[i].resetSteps();
    }
}

template<class ParamType, class PostProbCalculator, int BlockSize,
         class PostStepAction,class CustomDistribution, class Mover>
void EnsembleSampler<ParamType, PostProbCalculator, BlockSize, PostStepAction, CustomDistribution, Mover>::
    setSamplingMode(bool useSubSampling=false, int subSamplingInt=1, int burnIn=0)
{
    subSampling=useSubSampling;
    subSamplingInterval = subSamplingInt;
    markovChain.resetChainForSubSampling(burnIn, subSamplingInt);
    int tempSteps = (storedSteps - burnIn);
    storedSteps = (tempSteps/subSamplingInterval);
}

template<class ParamType, class PostProbCalculator, int BlockSize,
         class PostStepAction,class CustomDistribution, class Mover>
void EnsembleSampler<ParamType, PostProbCalculator, BlockSize, PostStepAction, CustomDistribution, Mover>::performStep(bool save)
{
    //first update all the red set from the black set
    for(int i=0; i<walkersPerSet; ++i)
    {
        //first get a proposed new point to move to
        ParamType scaling = moveProposer.getProposal(proposalPoint, numParams, walkerRedSet[i], walkerBlkSet, walkersPerSet, prng);
        walkerRedSet[i].proposePoint(proposalPoint, scaling, calc, prng, save);
    }
    //now update all the black set from the red set
    for(int i=0; i<walkersPerSet; ++i)
    {
        //first get a proposed new point to move to
        ParamType scaling = moveProposer.getProposal(proposalPoint, numParams, walkerBlkSet[i], walkerRedSet, walkersPerSet, prng);
        walkerBlkSet[i].proposePoint(proposalPoint, scaling, calc, prng, save);
    }
}

}
#endif  //MCMC_ENSEMBLESAMPLER_H
