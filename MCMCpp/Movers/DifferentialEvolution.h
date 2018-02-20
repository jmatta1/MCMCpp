/*!*****************************************************************************
********************************************************************************
**
** @copyright Copyright (C) 2018 James Till Matta
**
** This program is distributed in the hope that it will be useful,
** but WITHOUT ANY WARRANTY; without even the implied warranty of
** MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.
** 
********************************************************************************
*******************************************************************************/
#ifndef MCMCPP_MOVERS_DIFFERENTIALEVOLUTION_H
#define MCMCPP_MOVERS_DIFFERENTIALEVOLUTION_H
// includes for C system headers
// includes for C++ system headers
#include<cmath>
// includes from other libraries
// includes from MCMCpp..
#include"../Walker/Walker.h"
#include"../Utility/MultiSampler.h"
#include"../Utility/GwDistribution.h"
#include"../Utility/UserOjbectsTest.h"
#include"../Utility/Misc.h"

namespace MCMC
{
namespace Mover
{
/**
 * @class DifferentialEvolution
 * @ingroup Movers
 * @brief An object that calculates the next proposed step for a walker using the differential evolution algorithm
 * @author James Till Matta
 * 
 * \tparam ParamType The floating point type to be used for the chain, float, double, long double, etc.
 * \tparam Calculator The class that calculates the log posterior and whatever else a mover may need
 * 
 * A mover that uses the differential evolution MCMC technique presented in C.J.F. Ter Braak. Stat Comput (2006) 16:239–249
 */
template <class ParamType, class Calculator>
class DifferentialEvolution
{
    typedef Utility::GwDistribution<ParamType, 2, 1> DistType;
public:
    typedef Walker::Walker<ParamType> WalkType;
    static_assert(Utility::CheckCalcLogPostProb<Calculator, ParamType, ParamType*>(),
                  "DifferentialEvolution: The Calculator class does not have the necessary member function with signature:\n"
                  "  'ParamType calcLogPostProb(ParamType* paramSet)'");
    static_assert(std::is_copy_constructible<Calculator>::value, "The Calculator class needs to be copy constructible.");
    
    /*!
     * \brief Constructs a new stretch move object
     * \param numParams The number of parameters to work with
     * \param prngInit The seed for the random number generator
     * \param orig The original calculator class that will be copied to make the one stored internally
     */
    DifferentialEvolution(int numParams, long long prngInit, const Calculator& orig):
        paramCount(numParams), gamma(2.38/std::sqrt(2*paramCount)), prng(prngInit), calc(orig)
    {
        
        proposal = Utility::autoAlignedAlloc<ParamType>(paramCount);
    }
    
    ~DifferentialEvolution(){Utility::delAAA(proposal);}
    
    /*!
     * \brief Copy constructor
     * \param rhs Original DifferentialEvolution object to be copied
     */
    DifferentialEvolution(const DifferentialEvolution<ParamType, Calculator>& rhs):
    paramCount(rhs.paramCount), gamma(2.38/std::sqrt(2*paramCount)), calc(rhs.calc)
    {
        proposal = Utility::autoAlignedAlloc<ParamType>(paramCount);
    }
    
    /*!
     * \brief Deleted assignment operator
     */
    DifferentialEvolution<ParamType, Calculator>& operator=(const DifferentialEvolution<ParamType, Calculator>& rhs) = delete;
    
    /*!
     * \brief Allows the user to override the default choice for the scaling factor gamma (2.38 / sqrt(2*numDim))
     * \param newGamma The new value for gamma
     * 
     * According to the paper, for a multivariate normal target 2.38/sqrt(2*numDims) is about optimal, but
     * this is not necessarily the optimal choice for different targets
     */
    void overrideGamma(ParamType newGamma){ gamma = newGamma;}
    
    /*!
     * \brief Set a new b for the part of the proposal that needs UniformRand(-b, b)
     * \param newRand the new value for b
     * 
     * According to the paper, b can be quite small, infact in their numerical tests they used b = 10^-4 so that is my default
     */
    void overrideRandBounds(ParamType newRand){smallRandLowEdge = -newRand; smallRandWidth = 2*newRand;}
    
    /*!
     * \brief setPrngSeed Sets the seed and stream number of the underlying prng
     * \param seed The seed for the prng
     * \param stream The stream number for the prng
     */
    void setPrng(long long seed, long long stream){prng.setPrng(seed, stream);}
    
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
        //first choose the two points in walker set to calculate the difference to jump
        int ind1 = prng.getNonOffSetInt(numWalkers);
        int ind2;
        do{
            ind2 = prng.getNonOffSetInt(numWalkers);
        }while(ind2 == ind1);
        //now get the current walker's state
        const ParamType* currState = currWalker.getCurrState();
        //now get the states of the two selected walkers
        const ParamType* w1State = walkerSet[ind1].getCurrState();
        const ParamType* w2State = walkerSet[ind2].getCurrState();
        //now calculate the newly proposed point
        for(int i=0; i<paramCount; ++i)
        {
            proposal[i] = currState[i] +
                          gamma*(w1State[i] - w2State[i]) +
                          prng.getUniformRangeReal(smallRandLowEdge, smallRandWidth);
        }
        ParamType currAuxVal = calc.calcLogPostProb(proposal);
        if((currAuxVal - currWalker.getCurrAuxData()) > (prng.getNegExponentialReal()))
        {
            currWalker.jumpToNewPointSwap(proposal, currAuxVal, storePoint);
        }
        else
        {
            currWalker.stayAtCurrentPoint(storePoint);
        }
    }

private:
    ParamType* proposal = nullptr;
    int paramCount;
    ParamType gamma;
    ParamType smallRandWidth = 2.0e-4;
    ParamType smallRandLowEdge = -1.0e-4;
    Utility::MultiSampler<ParamType, DistType> prng;
    Calculator calc;
};
}
}
#endif  //MCMCPP_MOVERS_DIFFERENTIALEVOLUTION_H
