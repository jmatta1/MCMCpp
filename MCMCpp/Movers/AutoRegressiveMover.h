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
 * 
 * A fast, efficient Affine Invariant move algorithm functor that uses minimal resources and can yield good autocorrelation times
 */
template <class ParamType>
class AutoRegressiveMove
{
public:
    typedef Walker::Walker<ParamType> WalkType;
    
    /*!
     * \brief StretchMove Constructs a new stretch move object
     * \param numParams The number of parameters to work with
     * \param prngInit The seed for the random number generator
     * \param orig The original calculator class that will be copied to make the one stored internally
     */
    AutoRegressiveMove(int numParams, const ParamType* offsets, const ParamType* phiValues,
                       const ParamType* paramVariance): paramCount(numParams)
    {
        proposal = new ParamType[paramCount];
        phis = new ParamType[paramCount];
        offs = new ParamType[paramCount];
        prngStdDev = new ParamType[paramCount];
        for(int i=0; i<paramCount; ++i)
        {
            phis[i] = phiValues[i];
            offs[i] = offsets[i];
            prngStdDev[i] = std::sqrt(paramVariance[i])*std::sqrt(static_cast<ParamType>(1)-phis[i]*phis[i]);
        }
    }
    
    ~AutoRegressiveMove(){delete[] proposal;}
    
    /*!
     * \brief setPrngSeed Sets the seed of the underlying prng
     * \param seed The seed for the prng
     */
    void setPrngSeed(long long seed){prng.setPrngSeed(seed);}
    
    /*!
     * \brief getProposal Takes the curent walker, a set of walkers to draw a target from and calculates, assumes that
     * currWalker in not in the set of walkers to select from
     * \param currWalker A reference to the current walker that we are generating a proposal for
     * \param walkerSet A pointer to the set of walkers used to generate the proposal this is ignored for this mover
     * \param numWalkers The number of walkers in WalkerSet this is ignored for this mover
     * \param storePoint False if the point should not be written into the chain, accepted or not
     */
    void getProposal(WalkType& currWalker, WalkType* walkerSet, int numWalkers, bool storePoint)
    {
        ParamType* currState = currWalker.getCurrState();
        for(int i=0; i<paramCount; ++i)
        {
            proposal[i] = (offs[i] + (phis[i]*currState[i]) + (prngStdDev[i]*prng.getNormalReal()));
        }
        currWalker.jumpToNewPoint(proposal, static_cast<ParamType>(0), storePoint);
    }

    /*!
     * \brief getInitialPoints Generates an array of initial points for the walkers with the appropriate variance
     * \param initArray The array to hold the initial points
     * \param numWalkers The number of walkers to generate initial points for
     */
    void getInitialPoints(ParamType* initArray, int numWalkers)
    {
        for(int i=0; i<paramCount; ++i)
        {
            ParamType totalStdDev = prngStdDev[i]/std::sqrt(static_cast<ParamType>(1)-phis[i]*phis[i]);
            for(int j=0; j<numWalkers; ++j)
            {
                initArray[j*paramCount + i] = totalStdDev*prng.getNormalReal();
            }
        }
    }
    
private:
    ParamType* proposal = nullptr; ///<Holds the array that will tell the walker its new point
    ParamType* phis = nullptr; ///<Holds the array of recursion parameters
    ParamType* offs = nullptr; ///<Holds the array of offsets to the AR(1) model)
    ParamType* prngStdDev = nullptr; ///<Holds the standard deviation of each parameters needed normal distributed random number
    
    int paramCount; ///<Holds the total number of parameters
    Utility::MultiSampler<Utility::GwDistribution<ParamType, 2, 1> prng;  ///<The random number generator with an unused distribution transformer
};
}
}
#endif  //MCMC_WALKER_MOVERS_STRETCHMOVE_H
