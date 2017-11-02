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
#ifndef MCMC_UTILITY_GWDISTRIBUTION_H
#define MCMC_UTILITY_GWDISTRIBUTION_H
// includes for C system headers
// includes for C++ system headers
// includes from other libraries
// includes from MCMC

namespace MarkovChainMonteCarlo
{

namespace Utility
{

/**
 * @class GwDistribution
 * @ingroup Utility
 * @brief The class used to convert a uniform random fp into an fp drawn from 1/sqrt(z) with parameter Alpha
 * @author James Till Matta
 * 
 * \tparam ParamType The floating point type returned in the uniform and non-uniform floating point distributions
 * \tparam AlphaNum The integer in the numerator of the rational representation of Alpha (across which the distribution ranges 1/alpha to alpha)
 * \tparam AlphaDenom The integer in the denomerator of the rational representation of Alpha (across which the distribution ranges 1/alpha to alpha)
 * 
 * This functor class converts a floating point number drawn from the uniform
 * distribution between 0 and 1 to a floating point number draw from the
 * distribution proportional to 1/sqrt(z) on the interval AlphaDenom/AlphaNum to AlphaNum/AlphaDenom
 * i.e. 1/Alpha to Alpha
 */
template <class ParamType, int AlphaNum, int AlphaDenom>
class GwDistribution
{
public:
    //These will all be dumped into compile time constants
    constexpr static ParamType Alpha = (static_cast<ParamType>(AlphaNum)/static_cast<ParamType>(AlphaDenom));
    constexpr static ParamType SqrtAlpha = (std::sqrt(Alpha));
    constexpr static ParamType InvSqrtAlpha = (static_cast<ParamType>(1)/std::sqrt(Alpha));
    constexpr static ParamType Term1 = (SqrtAlpha - InvSqrtAlpha);

    inline ParamType operator()(ParamType in){ParamType temp = (Term1*in + InvSqrtAlpha); return (temp*temp);}
private:
};

}
}
#endif  //MCMC_UTILITY_GWDISTRIBUTION_H
