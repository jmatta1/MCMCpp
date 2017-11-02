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
#include"../Walker/Walker.h"
#include"../Utility/MultiSampler.h"
#include"../Utility/GwDistribution.h"
#include"../Utility/UserOjbectsTest.h"

namespace MarkovChainMonteCarlo
{
namespace Mover
{
/**
 * @class StretchMove
 * @ingroup Movers
 * @brief An object that calculates the next proposed step for a walker using the stretch move algorithm
 * @author James Till Matta
 * 
 * \tparam ParamType The floating point type to be used for the chain, float, double, long double, etc.
 * \tparam Calculator The class that calculates the log posterior and whatever else a mover may need
 * \tparam CustomDistribution The custom distribution used to draw samples for the various kinds of moves and other needed random numbers
 * 
 * A fast, efficient Affine Invariant move algorithm functor that uses minimal resources and can yield good autocorrelation times
 */
template <class ParamType, class Calculator, class CustomDistribution=Utility::GwDistribution<ParamType, 2, 1> >
class StretchMove
{
public:
    typedef Walker<ParamType> WalkType;
    static_assert(Utility::CheckCalcLogPostProb<Calculator, ParamType, ParamType*>(),
                  "StretchMove: The Calculator class does not have the necessary member function with signature:\n"
                  "  'ParamType calcLogPostProb(ParamType* paramSet)'");
    static_assert(Utility::CheckFunctor<CustomDistribution, ParamType, ParamType>(),
                  "The CustomDistribution class does not have the necessary member function with signature:\n"
                  "  'ParamType operator()(ParamType)'");
    static_assert(std::is_copy_constructible<Calculator>::value, "The Calculator class needs to be copy constructible.");
    static_assert(std::is_trivially_constructible<CustomDistribution>::value, "The CustomDistribution class needs to be trivially constructible.");
    
    /*!
     * \brief StretchMove Constructs a new stretch move object
     * \param numParams The number of parameters to work with
     * \param prngInit The seed for the random number generator
     * \param orig The original calculator class that will be copied to make the one stored internally
     */
    StretchMove(int numParams, long long prngInit, const Calculator& orig):
        paramCount(numParams), prng(prngInit), calc(orig)
    {proposal = new ParamType[paramCount];}
    
    ~StretchMove(){delete[] proposal;}
    
    StrethMove(const StretchMove<ParamType, Calculator, CustomDistribution>& rhs):
    paramCount(rhs.paramCount), calc(rhs.calc)
    {
        proposal = new ParamType[paramCount];
    }
    
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
        const ParamType* selectedState = walkerSet[prng.getNonOffSetInt(numWalkers)].getCurrState();
        const ParamType* currState = currWalker.getCurrState();
        ParamType scalingFactor = prng.getCustomSample();
        for(int i=0; i<numParams; ++i)
        {
            proposal[i] = (selectedState[i] + scalingFactor*(currState[i] - selectedState[i]));
        }
        //now calculate the log of the probability of performing the jump
        ParamType probScaling = (std::log(scalingFactor)*static_cast<ParamType>(numParams - 1));
        ParamType newProb = calc.calcLogPostProb(proposal);
        ParamType logProbDiff = (probScaling + newProb - currWalker.getCurrAuxData());
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
    ParamType* proposal = nullptr;
    int paramCount;
    Utility::MultiSampler<ParamType, CustomDistribution> prng;
    Calculator calc;
};
}
}
#endif  //MCMC_WALKER_MOVERS_STRETCHMOVE_H
