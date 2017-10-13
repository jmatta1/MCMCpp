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
// includes from other libraries
// includes from MCMC
#include"Chain/Chain.h"
#include"Walker/Movers/StretchMove.h"
#include"Walker/Movers/WalkMove.h"
#include"Walker/Walker.h"
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
 * @tparam LikelihoodCalculator The calculator of the log of the posterior probability function.
 * This can be calculated as the log of the likelihood function summed with the log of the prior-probability.
 * @tparam Mover The mover class, current options are StretchMove and WalkMove, StretchMove is
 * computationally much faster, but in some cases it can result in longer auto-correllation times,
 * WalkMove is much more computationally intensive, but it can result in shorter auto-correllation times in some cases
 * @tparam PostStepAction A functor that is called and given a pointer to the chain at the end of every step,
 * it need not be concurrent, it is called by a single thread, though that thread may be different for each call
 * @tparam PostBlockAction A functor that is called and given a pointer to the chain at the filling of every block
 * it need not be concurrent, it is called by a single thread, though that thread may be different for each call
 * @tparam CustomDistribution The class that will be used to set the get sample from custom
 * distribution function in the MultiSampler, currently only affects StretchMove
 * @tparam BlockSize Number of steps to store in a block of the linked list that stores the chain
 */
template<class ParamType,
         class LikelihoodCalculator,
         class Mover=Walker::StretchMove,
         class CustomDistribution=Utility::GwDistribution<ParamType, 2.0>,
         int BlockSize=1000,
         class PostStepAction=Utility::NoAction<ParamType, BlockSize>,
         class PostBlockAction=Utility::NoAction<ParamType, BlockSize>>
class EnsembleSampler
{
public:
    //define some useful typenames
    using ChainType=Chain::Chain<ParamType, BlockSize>;
    using MoverType=Mover<ParamType, BlockSize, CustomDistribution, LikelihoodCalculator>;
    using WalkerType=Walker::Walker<ParamType, BlockSize, MoverType, LikelihoodCalculator>;
    using PrngType=Utility::MultiSampler<ParamType, CustomDistribution>;
    using EnsembleSamplerType=EnsembleSampler<ParamType, LikelihoodCalculator, PostStepAction, PostBlockAction, Mover, CustomDistribution>;
    //perform static checks of the users classes to ensure that they have the needed member functions for their role
    static_assert(Utility::CheckPerformAction<PostStepAction, void, EnsembleSamplerType>(),
                  "The PostStepAction class does not have the necessary member function with signature:\n"
                  "  'void PerformAction(EnsembleSampler<ParamType, LikelihoodCalculator, PostStepAction,"
                  "                                      PostBlockAction, Mover, CustomDistribution>*)'");
    static_assert(Utility::CheckPerformAction<PostBlockAction, void, EnsembleSamplerType>(),
                  "The PostBlockAction class does not have the necessary member function with signature:\n"
                  "  'void PerformAction(EnsembleSampler<ParamType, LikelihoodCalculator, PostStepAction,\n"
                  "                                      PostBlockAction, Mover, CustomDistribution>*)'");
    static_assert(std::is_trivially_constructible<CustomDistribution>::value, "The custom distribution needs to be trivially constructible");
    /*!
     * \brief EnsembleSampler Constructs the ensemble sampler
     * \param runNumber RunNumber, used to seed the random number generator
     * \param numWalker The number of walkers to use
     * \param numParameter The number of parameters in the problem
     * \param maxChainSizeBytes The maximum size of the sample chain in bytes
     */
    EnsembleSampler(int runNumber, int numWalker, int numParameter, unsigned long long maxChainSizeBytes);
    
    /*!
     * @brief ~EnsembleSampler Delete the walker list then allow the rest to die normally
     */
    ~EnsembleSampler(){delete[] walkerList;}
    
    /*!
     * \brief setInitialWalkerPos Gives an initial position to every walker.
     * \param positions An array of floating point types with length numParameters*numWalker
     */
    void setInitialWalkerPos(ParamType* positions);
    
    /*!
     * \brief runMCMC Runs the Markov chain Monte Carlo for a set number of samples
     * \param numSamples The number of samples to store
     * 
     * In normal sampling mode, this function will run the ensemble for numSamples per walker
     * In subSampling mode, this function will run the ensemble numSamples*subSampingInterval and store numSamples
     */
    void runMCMC(int numSamples);
    
    /*!
     * \brief setSamplingMode Is used to change the sampling mode, either using subsampling, or taking all samples
     * \param useSubSampling Boolean for whether or not to use subsampling, default is false
     * \param subSamplingInt The interval for subsampling
     * \param burnIn The number of points at the beginning to discard for burnin
     */
    void setSamplingMode(bool useSubSampling=false, int subSampingInt=1, int burnIn=0);
    
    /*!
     * \brief getParamSetIttBegin Gets an iterator pointing to the beginning of the chain
     * that traverses the chain individual parameter set by individual parameter set. Incrementing
     * the iterator gets you the parameter set of the next walker, not paying attention to which parameter set
     * belongs to which step
     * \return A Parameter Set iterator pointed to the beginning of the chain.
     */
    Chain::ChainPsetIterator<ParamType, BlockSize> getParamSetIttBegin();
    /*!
     * \brief getParamSetIttEnd Gets a parameter set iterator pointing to the end of the chain
     * \return A parameter set iterator that points just past the end of the chain
     */
    Chain::ChainPsetIterator<ParamType, BlockSize> getParamSetIttEnd();
    
    /*!
     * \brief getStepIttBegin Gets a step iterator pointing to the chain step by step.
     * Incrementing the iterator takes you to the first walker of the next step
     * \return A Step iterator pointed to the beginning of the chain
     */
    Chain::ChainStepIterator<ParamType, BlockSize> getStepIttBegin();
    /*!
     * \brief getStepIttEnd Gets a step iterator that points to the end of chain
     * \return A Step iterator pointed to the end of the chain
     */
    Chain::ChainStepIterator<ParamType, BlockSize> getStepIttEnd();
    
    /*!
     * \brief getAutoCorrelationTimes Calculates the correlation time for a given parameter
     * \param paramIndex The index of the parameter [0, numParameter)
     * \param minAutoCorrTimes The minimum required number of autocorrelation times the algorithm needs to examine
     * \param step The increase in window size for each iteration of the algorithm
     * \param loWin The minimum window size
     * \param hiWin The maximum window size
     * \param fast If true, only use the first power of two samples to accellerate calculation using FFT
     * \return The autocorrelation time in samples for parameter # paramIndex
     */
    ParamType getAutoCorrelationTimes(int paramIndex, int minAutoCorrTimes=10, int step=1, int loWin=10, int hiWin=10000, bool fast=false);
private:
    ChainType markovChain;
    WalkerType* walkerList;
    MoverType moveProposer;
    PrngType prng;
    
    int numParams;
    int numWalkers;
    int subSampingInterval=1;
    
    bool isInitted = false;
    bool subSampling = false;
};

}
#endif  //MCMC_ENSEMBLESAMPLER_H
