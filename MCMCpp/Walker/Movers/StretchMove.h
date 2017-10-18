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
#include<cmath>
// includes from other libraries
// includes from MCMC
#include"../../Utility/MultiSampler.h"
#include"../Walker.h"

namespace MarkovChainMonteCarlo
{
namespace Walker
{
/**
 * @class StretchMove
 * @ingroup Walker
 * @ingroup Movers
 * @brief An object that calculates the next proposed step for a walker using the stretch move algorithm
 * @author James Till Matta
 * 
 * \tparam ParamType The floating point type to be used for the chain, float, double, long double, etc.
 * \tparam BlockSize The number of steps per walker that the block will hold
 * \tparam CustomDistribution The custom distribution used to draw samples for the various kinds of moves and other needed random numbers
 * \tparam PostProbCalculator The class that calculates the log posterior probability
 * 
 * A fast, efficient Affine Invariant move algorithm functor that uses minimal resources and can yield good autocorrelation times
 */
template <class ParamType, int BlockSize, class CustomDistribution, class PostProbCalculator>
class StretchMove
{
public:
    typedef Walker<ParamType, BlockSize, CustomDistribution, PostProbCalculator> WalkType;
    StretchMove(int numParams){} ///<useless constructor, exists to give the same constructor signature as WalkMove
    
    /*!
     * \brief getProposal Takes the curent walker, a set of walkers to draw a target from and calculates, assumes that
     * currWalker in not in the set of walkers to select from
     * \param proposal An array of ParamTypes that holds the proposed point this code will generate
     * \param numParams Number of parameters in the proposal array
     * \param currWalker A reference to the current walker that we are generating a proposal for
     * \param walkerSet A pointer to the set of walkers used to generate the proposal
     * \param numWalkers The number of walkers in WalkerSet
     * \param prng The pseudo random number generator to first select another walker for the calculation and to generate the scaling factor
     * \return The scaling factor for the post-prob ratio to be constructed for determining if the move should be taken
     * (z^(n-1) for this algorithm, but since we are working in logs it is (n-1)*z
     */
    ParamType getProposal(ParamType* proposal, int numParams, WalkType& currWalker, WalkType* walkerSet, int numWalkers, Utility::MultiSampler<ParamType, CustomDistribution>& prng)
    {
        WalkType& selWalker = walkerSet[prng.getNonOffSetInt(numWalkers)];
        ParamType scalingFactor = prng.getCustomSample();
        for(int i=0; i<numParams; ++i)
        {
            proposal[i] = (selWalker.currState[i] + scalingFactor*(currWalker.currState[i] - selWalker.currState[i]));
        }
        //now calculate the probability scaling factor
        return (std::log(scalingFactor)*static_cast<ParamType>((numParams - 1)));
    }

private:
};
}
}
#endif  //MCMC_WALKER_MOVERS_STRETCHMOVE_H
