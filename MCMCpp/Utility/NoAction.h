/*!*****************************************************************************
********************************************************************************
**
** @copyright Copyright (C) 2017-2018 James Till Matta
**
** This program is distributed in the hope that it will be useful,
** but WITHOUT ANY WARRANTY; without even the implied warranty of
** MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.
** 
********************************************************************************
*******************************************************************************/
#ifndef MCMCPP_UTILITY_NOACTION_H
#define MCMCPP_UTILITY_NOACTION_H
// includes for C system headers
// includes for C++ system headers
// includes from other libraries
// includes from MCMCpp
#include"../Chain/ChainStepIterator.h"

namespace MCMC
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
 */
template <class ParamType>
class NoAction
{
public:
    /*!
     * \brief PerformAction The perform action function, takes iterators to the begin and end of the chain and does *nothing*
     * \param startItt The iterator to the beginning of the chain
     * \param endItt The iterator to the end of the chain
     */
    void performAction(const Chain::ChainStepIterator<ParamType>& startItt,
                       const Chain::ChainStepIterator<ParamType>& endItt){}
private:
};

}
}
#endif  //MCMCPP_UTILITY_NOACTION_H
