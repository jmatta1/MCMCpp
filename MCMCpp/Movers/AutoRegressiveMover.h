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
#include<memory>
#include<iostream>
// includes from other libraries
// includes from MCMC
#include"../Walker/Walker.h"
#include"../Utility/ArrayDeleter.h"
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
 * @brief An object that calculates the next proposed step for a walker by using an AR1 model, useful for testing
 * @author James Till Matta
 * 
 * \tparam ParamType The floating point type to be used for the chain, float, double, long double, etc.
 * 
 * A mover that generates each next point using an autoregression model
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
    AutoRegressiveMove(int numParams, const ParamType* offsets, const ParamType* phiValues, const ParamType* paramVariance):
        paramCount(numParams), proposal(new ParamType[numParams]),
        phis(new ParamType[numParams], Utility::ArrayDeleter<ParamType>()),
        offs(new ParamType[numParams], Utility::ArrayDeleter<ParamType>()),
        prngStdDev(new ParamType[numParams], Utility::ArrayDeleter<ParamType>())
    {
        for(int i=0; i<paramCount; ++i)
        {
            phis.get()[i] = phiValues[i];
            offs.get()[i] = offsets[i];
            prngStdDev.get()[i] = std::sqrt(paramVariance[i])*std::sqrt(static_cast<ParamType>(1)-phiValues[i]*phiValues[i]);
        }
    }
    
    ~AutoRegressiveMove(){delete[] proposal;}
    
    AutoRegressiveMove(const AutoRegressiveMove<ParamType>& rhs):
        paramCount(rhs.paramCount), proposal(new ParamType[rhs.paramCount]),
        phis(rhs.phis), offs(rhs.offs), prngStdDev(rhs.prngStdDev){}
    
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
    void updateWalker(WalkType& currWalker, WalkType* walkerSet, int numWalkers, bool storePoint)
    {
        const ParamType* currState = currWalker.getCurrState();
        for(int i=0; i<paramCount; ++i)
        {
            proposal[i] = (offs.get()[i] + (phis.get()[i]*currState[i]) + (prngStdDev.get()[i]*prng.getNormalReal()));
        }
        currWalker.jumpToNewPoint(proposal, static_cast<ParamType>(0), storePoint);
    }

    /*!
     * \brief getInitialPoints Generates an array of initial points for the walkers with the appropriate variance
     * \param initArray The array to hold the initial points
     * \param numWalkers The number of walkers to generate initial points for
     */
    void getInitialPoints(ParamType* initArray, ParamType* auxArray, int numWalkers)
    {
        for(int i=0; i<paramCount; ++i)
        {
            auxArray[0] = static_cast<ParamType>(0);
            ParamType phi = phis.get()[i];
            ParamType totalStdDev = prngStdDev.get()[i]/std::sqrt(static_cast<ParamType>(1)-phi*phi);
            for(int j=0; j<numWalkers; ++j)
            {
                initArray[j*paramCount + i] = totalStdDev*prng.getNormalReal();
            }
        }
    }
    
private:
    int paramCount; ///<Holds the total number of parameters
    ParamType* proposal; ///Holds the new parameter set for the walker
    std::shared_ptr<ParamType> phis; ///<Holds the array of recursion parameters, doesn't need to be replicated
    std::shared_ptr<ParamType> offs; ///<Holds the array of offsets to the AR(1) model), doesn't need to be replicated
    std::shared_ptr<ParamType> prngStdDev; ///<Holds the standard deviation of each parameters needed normal distributed random number, doesn't need to be replicated
    Utility::MultiSampler<ParamType, Utility::GwDistribution<ParamType, 2, 1> > prng;  ///<The random number generator with an unused distribution transformer
};
}
}
#endif  //MCMC_WALKER_MOVERS_STRETCHMOVE_H
