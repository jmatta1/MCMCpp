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
#ifndef MCMCPP_UTILITY_GWDISTRIBUTION_H
#define MCMCPP_UTILITY_GWDISTRIBUTION_H
// includes for C system headers
// includes for C++ system headers
// includes from other libraries
// includes from MCMCpp

namespace MCMC
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
    //GCC seems to support some major mojo where it can use standard math functions in constexpr for c++11, not so much clang
#ifdef GCC_SUPPORTS_FANCYNESS
    constexpr static ParamType Alpha = (static_cast<ParamType>(AlphaNum)/static_cast<ParamType>(AlphaDenom));
    constexpr static ParamType SqrtAlpha = (std::sqrt(Alpha));
    constexpr static ParamType InvSqrtAlpha = (static_cast<ParamType>(1)/std::sqrt(Alpha));
    constexpr static ParamType Term1 = (SqrtAlpha - InvSqrtAlpha);
#else
    const ParamType Alpha = (static_cast<ParamType>(AlphaNum)/static_cast<ParamType>(AlphaDenom));
    const ParamType SqrtAlpha = (std::sqrt(Alpha));
    const ParamType InvSqrtAlpha = (static_cast<ParamType>(1)/std::sqrt(Alpha));
    const ParamType Term1 = (SqrtAlpha - InvSqrtAlpha);
#endif


    inline ParamType operator()(ParamType in){ParamType temp = (Term1*in + InvSqrtAlpha); return (temp*temp);}
private:
};

}
}
#endif  //MCMCPP_UTILITY_GWDISTRIBUTION_H
