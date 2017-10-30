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
#ifndef MCMC_WALKER_MOVERS_WALKMOVE_H
#define MCMC_WALKER_MOVERS_WALKMOVE_H
// includes for C system headers
// includes for C++ system headers
#include<cassert>
// includes from other libraries
// includes from MCMC
#include"../Utility/MultiSampler.h"
#include"../Utility/UserOjbectsTest.h"
#include"../Utility/GwDistribution.h"
#include"../Walker/Walker.h"

namespace MarkovChainMonteCarlo
{
namespace Mover
{
/**
 * @class WalkMove
 * @ingroup Movers
 * @brief An object that calculates the next proposed step for a walker using the walk move algorithm
 * @author James Till Matta
 * 
 * \tparam ParamType The floating point type to be used for the chain, float, double, long double, etc.
 * \tparam Calculator The class that calculates the log posterior and whatever else a mover may need
 * 
 * This mover is substantially more computationally expensive than StretchMove,
 * however in some circumstances it can give better autocorrellation times / chain properties
 */
template <class ParamType, class Calculator>
class WalkMove
{
public:
    typedef Walker::Walker<ParamType> WalkType;
    static_assert(Utility::CheckCalcLogPostProb<Calculator, ParamType, ParamType*>(),
                  "WalkMove: The Calculator class does not have the necessary member function with signature:\n"
                  "  'ParamType calcLogPostProb(ParamType* paramSet)'");
    static_assert(std::is_copy_constructible<Calculator>::value, "The Calculator class needs to be copy constructible.");
    
    /*!
     * \brief WalkMove Constructs the walk move object
     * \param numParams The number of parameters to work with
     * \param prngInit The seed for the random number generator
     * \param orig The original calculator class that will be copied to make the one stored internally
     * \param numSamples The number of samples to draw from the complementary ensemble
     */
    WalkMove(int numParams, long long prngInit, const Calculator& orig, int numSamples):
        numPoints(numSamples),ptCount(static_cast<ParamType>(numSamples)), paramCount(numParams), prng(prngInit), calc(orig)
    {
        walkerIndices = new int[numPoints];
        randoms = new ParamType[numPoints];
        proposal = new ParamType[paramCount];
    }
    
    ~WalkMove(){delete[] selectedWalkers; delete[] randoms;}
    
    /*!
     * \brief setPrngSeed Sets the seed of the underlying prng
     * \param seed The seed for the prng
     */
    void setPrngSeed(long long seed){prng.setPrngSeed(seed);}
    
    /*!
     * \brief getProposal Takes the curent walker, a set of walkers to draw a target from and calculates, assumes that
     * currWalker in not in the set of walkers to select from
     * \param currWalker A reference to the current walker that we are generating a proposal for
     * \param walkerSet A pointer to the set of walkers used to generate the proposal
     * \param numWalkers The number of walkers in WalkerSet
     * \param storePoint False if the point should not be written into the chain, accepted or not
     */
    void updateWalker(WalkType& currWalker, WalkType* walkerSet, int numWalkers, bool storePoint)
    {
        selectWalkers(numWalkers);
        calculateProposal(currWalker, walkerSet);
        ParamType newProb = calc.calcLogPostProb(proposal);
        ParamType logProbDiff = (newProb - currWalker.getCurrAuxData());
        if(prng.getNegExponentialReal() < logProbDiff)
        {
            currWalker.jumpToNewPoint(proposal, newProb, storePoint);
        }
        else
        {
            currWalker.stayAtCurrentPoint(storePoint);
        }
    }

private:
    /*!
     * \brief selectWalkers chooses a set of walkers to calculate the point proposal
     * \param numWalkers The number of walkers available to choose from
     * 
     * This algorithm is drawn directly from Knuth's semi-numerical algorithms.
     * There is a much more complicated, but faster algorithm available in the article
     * J. Vitter. "An Efficient Algorithm for Sequential Random Sampling"
     * ACM Transactions on Mathematical Software, 13(1), March 1987, 58-67.
     * but this one is good enough for now
     */
    void selectWalkers(int numWalkers)
    {
        assert(numPoints <= numWalkers);//probably unnecessary, but safety first
        int numSelected = 0;
        int totalInputExamined = 0;
        ParamType randUniform;
        while(numSelected < numPoints)
        {
            randUniform = prng.getUniformReal();
            if( ((numWalkers-totalInputExamined)*randUniform) >= (numWalkers-numSelected))
            {
                ++totalInputExamined;
            }
            else
            {
                walkerIndices[numSelected] = totalInputExamined;
                ++totalInputExamined;
                ++numSelected;
            }
        }
    }
    
    /*!
     * \brief calculateProposal Takes the set of selected walkers and generates the proposal point
     */
    void calculateProposal(WalkType& currWalker, WalkType* walkerSet)
    {
        // do the first parameter explicitly to store the random values
        // we do this rather than simply loop to generate the random values because this way we do numPoints*paramCount iterations
        // the other way requires (paramCount+1)*numPoints iterations, its probably small, but every little bit helps
        intermediates[0] = 0.0;
        intermediates[1] = 0.0;
        intermediates[2] = 0.0;
        for(int j=0; j<numPoints; ++j)
        {
            randoms[j] = prng.getNormalReal();
            ParamType value = walkerSet[walkerIndices[j]].getCurrState()[0];
            intermediates[0] += randoms[j]*value;
            intermediates[1] += randoms[j];
            intermediates[2] += value;
        }
        proposal[0] = currWalker.getCurrState()[0] + (intermediates[0] - (intermediates[1]*(intermediates[2]/ptCount)));
        //loop across the remaining parameters, using the pre-stored random numbers
        for(int i=1; i<paramCount; ++i)
        {
            intermediates[0] = 0.0;
            intermediates[1] = 0.0;
            intermediates[2] = 0.0;
            for(int j=0; j<numPoints; ++j)
            {
                ParamType value = walkerSet[walkerIndices[j]].getCurrState()[i];
                intermediates[0] += randoms[j]*value;
                intermediates[1] += randoms[j];
                intermediates[2] += value;
            }
            proposal[i] = (currWalker.getCurrState()[i] + (intermediates[0] - (intermediates[1]*(intermediates[2]/ptCount))));
        }
    }
    
    int numPoints; ///<Number of points to sample from to generate the point proposal
    ParamType ptCount; ///<Number of points to sample from to generate the point proposal stored as a double
    ParamType intermediates[3]; ///<Intermediate values for use in calculating the point proposal
    ParamType* randoms; ///<Storage for the random numbers selected in calculating the first parameter in the point proposal
    int* walkerIndices; ///<Storage for the indices of the randomly selected walkers
    
    ParamType* proposal = nullptr;
    int paramCount;
    Utility::MultiSampler<ParamType, Utility::GwDistribution<ParamType, 2, 1>> prng;
    Calculator calc;
};
}
}
#endif  //MCMC_WALKER_MOVERS_WALKMOVE_H
