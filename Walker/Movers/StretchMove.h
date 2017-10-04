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
 * @class StretchMove
 * @ingroup Walker
 * @ingroup Movers
 * @brief A functor object that calculates the next proposed step for a walker using the stretch move algorithm
 * @author James Till Matta
 * 
 * \tparam ParamType The floating point type to be used for the chain, float, double, long double, etc.
 * \tparam BlockSize The number of steps per walker that the block will hold
 * \tparam CustomDistribution The custom distribution used to draw samples for the various kinds of moves and other needed random numbers
 * \tparam LikelihoodCalculator The class that calculates the log likelihood
 * 
 */
template <class ParamType, int BlockSize, class CustomDistribution, class LikelihoodCalculator>
class StretchMove
{
public:
    typedef Walker<ParamType, BlockSize, CustomDistribution, LikelihoodCalculator> WalkType;
    
    /*!
     * \brief getProposal Takes the curent walker, a set of walkers to draw a target from and calculates
     * \param proposal An array of ParamTypes that holds the proposed point this code will generate
     * \param numParams Number of parameters in the proposal array
     * \param currWalker A reference to the current walker that we are generating a proposal for
     * \param walkerSet A pointer to the set of walkers used to generate the proposal
     * \param numWalkers The number of walkers in WalkerSet
     * \param prng The pseudo random number generator to first select another walker for the calculation and to generate the scaling factor
     * \return The scaling factor for the likelihood ratio to be constructed for determining if the move should be taken
     */
    ParamType getProposal(ParamType* proposal, int numParams, WalkType& currWalker, WalkType* walkerSet, int numWalkers, Utility::MultiSampler<CustomDistribution>& prng)
    {
        WalkType& selWalker = walkerSet[prng.getNonOffSetInt(numWalkers);]
        ParamType scalingFactor = prng.getCustomSample();
        for(int i=0; i<numParams; ++i)
        {
            proposal[i] = (selWalker.currState[i] + scalingFactor*(currWalker.currState[i] - selWalker.currState[i]));
        }
        //now do exponentiation by squaring on the scaling factor to return it
        int exponent = (numParams - 1);
        ParamType retValue = 1.0;
        while(exponent > 0)
        {
            if(exponent&0x00000001) retValue *= scalingFactor;
            exponent >>= 1;
            scalingFactor *= scalingFactor;
        }
        return retValue;
    }

private:
};
}
}
#endif  //MCMC_WALKER_MOVERS_STRETCHMOVE_H
