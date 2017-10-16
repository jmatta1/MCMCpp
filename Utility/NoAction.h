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
#ifndef MCMC_UTILITY_NOACTION_H
#define MCMC_UTILITY_NOACTION_H
// includes for C system headers
// includes for C++ system headers
// includes from other libraries
// includes from MCMC
#include"../Chain/ChainStepIterator.h"

namespace MarkovChainMonteCarlo
{

namespace Utility
{

/**
 * @class NoAction
 * @ingroup Utility
 * @brief This class serves the purpose of a end of * action class, that should, upon compilation, convert into an inlined nothing
 * @author James Till Matta
 * 
 * \tparam ParamType The floating point type returned in the uniform and non-uniform floating point distributions
 * \tparam BlockSize The number of steps stored in each block
 */
template <class ParamType, int BlockSize>
class NoAction
{
public:
    /*!
     * \brief PerformAction The perform action function, takes iterators to the begin and end of the chain and does *nothing*
     * \param startItt The iterator to the beginning of the chain
     * \param endItt The iterator to the end of the chain
     */
    void performAction(Chain::ChainStepIterator<ParamType, BlockSize>& startItt,
                       Chain::ChainStepIterator<ParamType, BlockSize>& endItt){}
private:
};

}
}
#endif  //MCMC_UTILITY_NOACTION_H
