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
#ifndef MCMC_ENSEMBLESAMPLER_H
#define MCMC_ENSEMBLESAMPLER_H
// includes for C system headers
// includes for C++ system headers
// includes from other libraries
// includes from MCMC
#include"Chain/Chain.h"
#include"Walker/Movers/StretchMove.h"
#include"Walker/Movers/WalkMove.h"
#include"Walker/Walker.h"
#include"Utility/GwDistribution.h"
#include"Utility/NoStepEndAction.h"

namespace MarkovChainMonteCarlo
{
/*!
 * @class EnsembleSampler
 * @ingroup Primary
 * @brief The user interface to perform a sequential MCMC sampling
 * @author James Till Matta
 * 
 * @tparam ParamType The floating point type for which this calculation will be carried out
 * @tparam LikelihoodCalculator The calculator of the log of the posterior probability function.
 * This can be calculated as the log of the likelihood function summed with the log of the prior-probability.
 * @tparam Mover The mover class, current options are StretchMove and WalkMove, StretchMove is
 * computationally much faster, but in some cases it can result in longer auto-correllation times,
 * WalkMove is much more computationally intensive, but it can result in shorter auto-correllation times in some cases
 * @tparam PostStepAction A functor that is called and given a pointer to the chain at the end of every step,
 * it need not be concurrent, it is called by a single thread, though that thread may be different for each call
 * @tparam PostBlockAction A functor that is called and given a pointer to the chain at the filling of every block
 * it need not be concurrent, it is called by a single thread, though that thread may be different for each call
 * @tparam CustomDistribution The class that will be used to set the get sample from custom
 * distribution function in the MultiSampler, currently only affects StretchMove
 * @tparam BlockSize Number of steps to store in a block of the linked list that stores the chain
 */

template<class ParamType,
         class LikelihoodCalculator,
         class PostStepAction=Utility::NoAction,
         class PostBlockAction=Utility::NoAction,
         class Mover=Walker::StretchMove,
         class CustomDistribution=Utility::GwDistribution<ParamType, 2.0>,
         int BlockSize=1000>
class EnsembleSampler
{
public:
    EnsembleSampler(int numWalker, int numParameter, unsigned long long maxChainSize);
    ~EnsembleSampler();
    
    void setInitialWalkerPos(ParamType* positions);
    
    void runMCMC(int numSamples);
    
    Chain::ChainPsetIterator<ParamType, BlockSize> getParamSetIttBegin();
    Chain::ChainPsetIterator<ParamType, BlockSize> getParamSetIttEnd();
    
    Chain::ChainStepIterator<ParamType, BlockSize> getStepIttBegin();
    Chain::ChainStepIterator<ParamType, BlockSize> getStepIttEnd();
    
    double getParameterAutoCorrelationTime(int paramIndex, int minAutoCorrTimes=10, int step=1, int loWin=10, int hiWin=10000, bool fast=false);
private:

};

}
#endif  //MCMC_ENSEMBLESAMPLER_H
