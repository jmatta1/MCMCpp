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
#include"../Chain/Chain.h"
#include"../Utility/MultiSampler.h"


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
 * \tparam PostProbCalculator The class that calculates the log post probability
 * 
 */
template <class ParamType, int BlockSize>
class Walker
{
public:
    /*!
     * \brief Walker Constructs a default instance of Walker with all values set to nullptr or zero
     */
    Walker(){}
    ~Walker(){if(currState != nullptr) delete[] currState;}
    
    /*!
     * \brief init Initializes the walker in question so that it can be used
     * \param chain A pointer to the chain that will store parameter sets
     * \param walkerIndex The index of the current walker
     * \param numParam The number of parameters in the spaces to be walked
     * \param numCell The number of cells to be stored per walker in the chain
     */
    void init(Chain::Chain<ParamType, BlockSize>* chain, int walkerIndex, int numParam)
    {markovChain = chain; walkerNumber = walkerIndex; numParams = numParam; currState = new ParamType[numParam];}
    
    /*!
     * \brief setFirstPoint Initializes the walker with its very first point (which is counted as an accepted step (and a normal step))
     * \param init The parameter set for the initial value
     * \param calc The log post prob calculator, see Walker::proposePoint for why it is passed and not stored
     */
    void setFirstPoint(ParamType* init, const ParamType& auxVal);
    
    /*!
     * \brief jumpToNewPoint Tells the walker there has been a new point accepted and that it should jump to it
     * \param newPos The new point for the walker
     * \param auxVal the value of the auxilliary data to be stored with the new position
     */
    void jumpToNewPoint(ParamType* newPos, const ParamType& auxVal);
    
    /*!
     * \brief jumpToNewPoint Tells the walker that the step to this point has been rejected and that it should stay at the current point
     */
    void stayAtCurrentPoint(){++totalSteps; markovChain->storeWalker(walkerNumber, currState);}
    
    /*!
     * \brief getTotalSteps Gets the total number of steps taken (initial position not counted)
     * \return The total number of steps taken
     */
    int getTotalSteps(){return totalSteps;}
    
    /*!
     * \brief getAcceptedProposals Gets the number of times a new proposal was taken (initial position not counted)
     * \return Gets the number of accepted proposals
     */
    int getAcceptedProposals(){return acceptedSteps;}
    
    /*!
     * \brief resetSteps Resets the accepted and total steps to 0
     */
    void resetSteps(){acceptedSteps = 0; totalSteps = 0;}
    
    /*!
     * \brief getCurrState returns the current point that the walker is at
     * \return array of size numParams containing the walker's current position
     */
    const ParamType* getCurrState(){return currState;}
    
    /*!
     * \brief getCurrAuxData returns the auxilliary data that is stored with the walker
     * \return the auxilliary data value
     */
    const ParamType getCurrAuxData(){return auxData;}
    
private:
    Chain::Chain<ParamType, BlockSize>* markovChain = nullptr; ///<Holds a reference to the chain that stores the points of the walker
    ParamType* currState = nullptr; ///<Holds the current position, should have at least one extra cell to hold the post prob for that position when written to the chain
    ParamType auxData = 0.0; ///< holds the auxiliary data the walker needs for a given position, to decide to make jumps or not
    int walkerNumber = 0; ///<Holds the integer index of the walker, unique between walkers
    int numParams = 0; ///<holds the number of parameters for the post prob function
    int acceptedSteps = 0; ///< holds the number of times that a new parameter set was accepted
    int totalSteps = 0; ///< holds the number of times that a new parameter set was proposed
};

template <class ParamType, int BlockSize>
void Walker<ParamType, BlockSize>::setFirstPoint(ParamType* init, const ParamType& auxVal)
{
    for(int i=0; i<numParams; ++i)
    {
        currState[i] = init[i];
    }
    auxData = auxVal;
    ++acceptedSteps;
    ++totalSteps;
    markovChain->storeWalker(walkerNumber, currState);
}

template <class ParamType, int BlockSize>
void Walker<ParamType, BlockSize>::jumpToNewPoint(ParamType* newPos, const ParamType& auxVal)
{
    for(int i=0; i<numParams; ++i)
    {
        currState[i] = newPos[i];
    }
    auxData = auxVal;
    ++acceptedSteps;
    ++totalSteps;
    markovChain->storeWalker(walkerNumber, currState);
}

}
}
#endif  //MCMC_WALKER_WALKER_H
