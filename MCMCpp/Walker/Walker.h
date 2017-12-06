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
#include<stdlib.h>//needed for aligned allocation, which appears in C11 but does not appear in C++ until C++17
// includes for C++ system headers
#include<algorithm>
// includes from other libraries
// includes from MCMC
#include"../Chain/Chain.h"
#include"../Utility/MultiSampler.h"
#include"Utility/Misc.h"


namespace MCMC
{

namespace Mover
{

}

namespace Walker
{
/**
 * @class Walker
 * @ingroup Walker
 * @brief The object that represents the state and actions of a single walker for the sequential algorithm
 * @author James Till Matta
 * 
 * \tparam ParamType The floating point type to be used for the chain, float, double, long double, etc.
 * 
 * **WARNING** If using the swap procedure instead of copy (for speed), movers need to allocate their proposal arrays
 * using std::aligned_alloc in cstdlib instead of using operator new
 */
template <class ParamType>
class Walker
{
public:
    /*!
     * \brief Walker Constructs a default instance of Walker with all values set to nullptr or zero
     */
    Walker(){}
    ~Walker(){if(currState != nullptr) free(currState);}
    
    /*!
     * \brief Delete copy constructor
     */
    Walker(const Walker<ParamType>& rhs) = delete;
    
    /*!
     * \brief Delete assignment operator
     */
    Walker<ParamType>& operator=(const Walker<ParamType>& rhs) = delete;
    
    /*!
     * \brief init Initializes the walker in question so that it can be used
     * \param chain A pointer to the chain that will store parameter sets
     * \param walkerIndex The index of the current walker
     * \param numParam The number of parameters in the spaces to be walked
     * \param numCell The number of cells to be stored per walker in the chain
     */
    void init(Chain::Chain<ParamType>* chain, int walkerIndex, int numParam)
    {
        markovChain = chain;
        walkerNumber = walkerIndex; 
        numParams = numParam;
        size_t allocSize = (sizeof(ParamType)*numParam);
        if(allocSize%Utility::AlignmentLength) //if allocSize is not an integral multiple of AlignementLength
        {
            allocSize = (((allocSize/Utility::AlignmentLength)+1)*Utility::AlignmentLength);
        }
        currState = reinterpret_cast<ParamType*>(aligned_alloc(Utility::AlignmentLength,allocSize));
    }
    
    /*!
     * \brief setFirstPoint Initializes the walker with its very first point (which is counted as an accepted step (and a normal step))
     * \param init The parameter set for the initial value
     * \param calc The log post prob calculator, see Walker::proposePoint for why it is passed and not stored
     */
    void setFirstPoint(ParamType* init, const ParamType& auxVal, bool storePoint=true)
    {
        jumpToNewPoint(init, auxVal, storePoint);
    }
    
    /*!
     * \brief storeCurrentPoint Stores the walker's current position into the chain, no questions asked
     */
    void storeCurrentPoint()
    {
        markovChain->storeWalker(walkerNumber, currState);
    }
    
    /*!
     * \brief jumpToNewPoint Tells the walker there has been a new point accepted and that it should jump to it, this version copies values
     * \param newPos The new point for the walker
     * \param auxVal the value of the auxilliary data to be stored with the new position
     * \param storePoint If false, do not store the new point in the Markov chain
     */
    void jumpToNewPoint(ParamType* newPos, const ParamType& auxVal, bool storePoint=true);
    
    /*!
     * \brief jumpToNewPoint Tells the walker there has been a new point accepted and that it should jump to it, this version swaps pointers with the mover
     * \param newPos The new point for the walker
     * \param auxVal the value of the auxilliary data to be stored with the new position
     * \param storePoint If false, do not store the new point in the Markov chain
     * 
     * \note If you use this, your newPos pointer needs to have been allocated using the c++ new[] operator, not calloc, malloc, or aligned_alloc
     */
    void jumpToNewPointSwap(ParamType*& newPos, const ParamType& auxVal, bool storePoint=true);
    
    /*!
     * \brief jumpToNewPoint Tells the walker that the step to this point has been rejected and that it should stay at the current point
     * \param storePoint If false, do not store the new point in the Markov chain
     */
    void stayAtCurrentPoint(bool storePoint=true){++diffSteps; if(storePoint) markovChain->storeWalker(walkerNumber, currState);}
    
    /*!
     * \brief getTotalSteps Gets the total number of steps taken (initial position not counted)
     * \return The total number of steps taken
     */
    int getTotalSteps(){return (acceptedSteps+diffSteps);}
    
    /*!
     * \brief getAcceptedProposals Gets the number of times a new proposal was taken (initial position not counted)
     * \return Gets the number of accepted proposals
     */
    int getAcceptedProposals(){return acceptedSteps;}
    
    /*!
     * \brief resetSteps Resets the accepted and total steps to 0
     */
    void resetSteps(){acceptedSteps = 0; diffSteps = 0;}
    
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
    
    int getWalkerNum(){return walkerNumber;}
private:
    Chain::Chain<ParamType>* markovChain = nullptr; ///<Holds a reference to the chain that stores the points of the walker
    ParamType* currState = nullptr; ///<Holds the current position, should have at least one extra cell to hold the post prob for that position when written to the chain
    ParamType auxData = 0.0; ///< holds the auxiliary data the walker needs for a given position, to decide to make jumps or not
    int walkerNumber = 0; ///<Holds the integer index of the walker, unique between walkers
    int numParams = 0; ///<holds the number of parameters for the post prob function
    int acceptedSteps = 0; ///< holds the number of times that a new parameter set was accepted
    int diffSteps = 0; ///< holds the number of times that a new parameter set was proposed and not accepted
};

template <class ParamType>
void Walker<ParamType>::jumpToNewPoint(ParamType* newPos, const ParamType& auxVal, bool storePoint)
{
    std::copy(newPos, (newPos+numParams), currState);
    auxData = auxVal;
    ++acceptedSteps;
    if(storePoint) markovChain->storeWalker(walkerNumber, currState);
}

template <class ParamType>
void Walker<ParamType>::jumpToNewPointSwap(ParamType* &newPos, const ParamType& auxVal, bool storePoint)
{
    std::swap(currState, newPos);
    auxData = auxVal;
    ++acceptedSteps;
    if(storePoint) markovChain->storeWalker(walkerNumber, currState);
}

}
}
#endif  //MCMC_WALKER_WALKER_H
