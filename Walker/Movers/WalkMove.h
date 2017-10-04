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
#ifndef MCMC_WALKER_MOVERS_STRETCHMOVE_H
#define MCMC_WALKER_MOVERS_STRETCHMOVE_H
// includes for C system headers
// includes for C++ system headers
// includes from other libraries
// includes from MCMC
#include"Reconstruction/MCMC/Utility/MultiSampler.h"
#include"../Walker.h"

namespace MarkovChainMonteCarlo
{
namespace Walker
{
/**
 * @class WalkMove
 * @ingroup Walker
 * @ingroup Movers
 * @brief An object that calculates the next proposed step for a walker using the walk move algorithm
 * @author James Till Matta
 * 
 * \tparam ParamType The floating point type to be used for the chain, float, double, long double, etc.
 * \tparam BlockSize The number of steps per walker that the block will hold
 * \tparam CustomDistribution The custom distribution used to draw samples for the various kinds of moves and other needed random numbers
 * \tparam LikelihoodCalculator The class that calculates the log likelihood
 * 
 * This mover is substantially more computationally expensive than StretchMove, however in some circumstances it can give better autocorrellation times
 */
template <class ParamType, int BlockSize, class CustomDistribution, class LikelihoodCalculator>
class WalkMove
{
public:
    WalkMove(int numPts):numPoints(numPts),ptCount(static_cast<ParamType>(numPts)){selectedWalkers = new int[numPoints]; randoms = new ParamType[numPoints];}
    ~WalkMove(){delete[] selectedWalkers; delete[] randoms;}
    typedef Walker<ParamType, BlockSize, CustomDistribution, LikelihoodCalculator> WalkType;
    
    /*!
     * \brief getProposal Takes the curent walker, a set of walkers to draw a target from and calculates
     * \param proposal An array of ParamTypes that holds the proposed point this code will generate
     * \param numParams Number of parameters in the proposal array
     * \param currWalker A reference to the current walker that we are generating a proposal for
     * \param walkerSet A pointer to the set of walkers used to generate the proposal
     * \param numWalkers The number of walkers in WalkerSet
     * \param prng The pseudo random number generator to first select another walker for the calculation and to generate the scaling factor
     * \return The scaling factor for the likelihood ratio to be constructed for determining if the move should be taken (1.0 for this algorithm)
     */
    ParamType getProposal(ParamType* proposal, int numParams, WalkType& currWalker, WalkType* walkerSet, int numWalkers, Utility::MultiSampler<ParamType, CustomDistribution>& prng)
    {
        selectWalkers(numWalkers, prng);
        // do the first parameter explicitly to store the random values
        intermediates[0] = 0.0;
        intermediates[1] = 0.0;
        intermediates[2] = 0.0;
        for(int j=0; j<numPoints; ++j)
        {
            randoms[j] = prng.getNormalReal();
            ParamType value = walkerSet[walkerIndices[j]].currState[0];
            intermediates[0] += randoms[j]*value;
            intermediates[1] += randoms[j];
            intermediates[2] += value;
        }
        proposal[0] = currWalker.currState[0] + (intermediates[0] - (intermediates[1]*(intermediates[2]/ptCount)));
        //loop across the remaining parameters, using the pre-stored random numbers
        for(int i=1; i<numParams; ++i)
        {
            intermediates[0] = 0.0;
            intermediates[1] = 0.0;
            intermediates[2] = 0.0;
            for(int j=0; j<numPoints; ++j)
            {
                ParamType value = walkerSet[walkerIndices[j]].currState[i];
                intermediates[0] += randoms[j]*value;
                intermediates[1] += randoms[j];
                intermediates[2] += value;
            }
            proposal[i] = (currWalker.currState[i] + (intermediates[0] - (intermediates[1]*(intermediates[2]/ptCount))));
        }
        return 1.0;
    }

private:
    void selectWalkers(int numWalkers, Utility::MultiSampler<CustomDistribution>& prng)
    {
        int i=0;
        while(i < numWalkers)
        {
            walkerIndices[i] = prng.getNonOffSetInt(numWalkers);
            bool noRepeats = true;
            for(int j = 0; j<i; ++j)
            {
                if(walkerIndices[i] == walkerIndices[j])
                {
                    noRepeats = false;
                }
            }
            if(noRepeats) ++i;
        }
    }
    
    int numPoints;
    ParamType ptCount;
    ParamType intermediates[3];
    ParamType* randoms;
    int* walkerIndices;
};
}
}
#endif  //MCMC_WALKER_MOVERS_STRETCHMOVE_H
