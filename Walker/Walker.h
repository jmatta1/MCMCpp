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
#ifndef MCMC_WALKER_WALKER_H
#define MCMC_WALKER_WALKER_H
// includes for C system headers
// includes for C++ system headers
// includes from other libraries
// includes from MCMC
#include"Reconstruction/MCMC/Chain/Chain.h"
#include"Reconstruction/MCMC/Utility/MultiSampler.h"


namespace MarkovChainMonteCarlo
{
namespace Walker
{
/**
 * @class Walker
 * @ingroup Walker
 * @brief The object that represents the state and actions of a single walker for the sequential algorithm
 * @author James Till Matta
 * 
 * \tparam ParamType The floating point type to be used for the chain, float, double, long double, etc.
 * \tparam BlockSize The number of steps per walker that the block will hold
 * \tparam CustomDistribution The custom distribution used to draw samples for the various kinds of moves and other needed random numbers
 * \tparam LikelihoodCalculator The class that calculates the log likelihood
 * 
 */
template <class ParamType, int BlockSize, class CustomDistribution, class LikelihoodCalculator>
class Walker
{
public:
    Walker(Chain::Chain<ParamType, BlockSize>* chain, int walkerIndex, int numParam, int numCell):
        markovChain(chain), walkerNumber(walkerIndex), numParams(numParam), numCells(numCell)
    {    currState = new ParamType[numCells];}
    ~Walker(){delete[] currState;}
    
    /*!
     * \brief setFirstPoint Initializes the walker with its very first point (which is counted as an accepted step (and a normal step))
     * \param init The parameter set for the initial value
     * \param calc The logLikelihood calculator, see Walker::proposePoint for why it is passed and not stored
     */
    void setFirstPoint(ParamType* init, LikelihoodCalculator& calc);
    
    /*!
     * \brief proposePoint Stores the walker's current point in the MarkovChainMonteCarlo::Chain::Chain object, checks the proposal point, and moves there if accepted
     * \param newPos The point that is proposed as a next point in the walker
     * \param ratioScale The scaling factor of the ratio for acceptance (Z^(n-1) for stretch move, 1.0 for walk move)
     * \param calc The likelihood calculator, this is passed instead of stored so that a thread can pass its local copy
     * \param prng The pseudorandom number generator sampler, passed instead of stored so that a thread can pass its local copy
     */
    void proposePoint(ParamType* newPos, const ParamType& ratioScale, LikelihoodCalculator& calc, Utility::MultiSampler<CustomDistribution>& prng);
    
    /*!
     * \brief saveFinalPoint Places the point currently stored on the walker into the chain without testing a new one.
     */
    void saveFinalPoint(){currState[numParams] = currLikelihood; markovChain->storeWalker(walkerNumber, currState);}
    
    /*!
     * \brief getTotalSteps Gets the total number of steps taken (Number of proposals + 1 because setting the walker initial state counts as a step)
     * \return The total number of steps taken
     */
    int getTotalSteps(){return totalSteps;}
    
    /*!
     * \brief getAcceptedProposals Gets the number of times a new proposal was taken (+1 because the initial state of the walker is counted as an accepted proposal)
     * \return Gets the number of accepted proposals
     */
    int getAcceptedProposals(){return acceptedSteps;}
    
private:
    Chain::Chain<ParamType, BlockSize>* markovChain; ///<Holds a reference to the chain that stores the points of the walker
    ParamType* currState = nullptr; ///<Holds the current position, should have at least one extra cell to hold the likelihood for that position when written to the chain
    ParamType currLikelihood; ///< holds the current likelihood (transfered into currState prior to dump to chain)
    int walkerNumber; ///<Holds the integer index of the walker
    int numParams; ///<holds the number of parameters for the likelihood function
    int numCells; ///<holds the number of cells in the currState array
    int acceptedSteps = 0; ///< holds the number of times that a new parameter set was accepted
    int totalSteps = 0; ///< holds the number of times that a new parameter set was proposed
};

template <class ParamType, int BlockSize, class CustomDistribution, class LikelihoodCalculator>
void Walker<ParamType, BlockSize, CustomDistribution, LikelihoodCalculator>::setFirstPoint(ParamType* init, LikelihoodCalculator& calc)
{
    for(int i=0; i<numParams; ++i)
    {
        currState[i] = init[i];
    }
    currLikelihood = calc.CalculateLogLikelihood(init);
    ++acceptedSteps;
    ++totalSteps;
}

template <class ParamType, int BlockSize, class CustomDistribution, class LikelihoodCalculator>
void Walker<ParamType, BlockSize, CustomDistribution, LikelihoodCalculator>::proposePoint(ParamType* newPos, const MarkovChainMonteCarlo::ParamType& ratioScale, LikelihoodCalculator& calc, Utility::MultiSampler<CustomDistribution>& prng)
{
    currState[numParams] = currLikelihood;
    markovChain->storeWalker(walkerNumber, currState);
    ParamType newLikelihood = calc.CalculateLogLikelihood(newPos);
    ParamType ratio = ratioScale*(newLikelihood/currLikelihood);
    if(prng.getUniformReal() < ratio)
    {
        for(int i=0; i<numParams; ++i)
        {
            currState[i] = newPos[i];
        }
        currLikelihood = newLikelihood;
        ++acceptedSteps;
    }
    ++totalSteps;
}

}
}
#endif  //MCMC_WALKER_WALKER_H
