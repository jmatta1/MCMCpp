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
#ifndef MCMC_UTILITY_MULTISAMPLER_H
#define MCMC_UTILITY_MULTISAMPLER_H
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
 * @brief The class used to convert a uniform random fp into an fp drawn from 1/sqrt(z)
 * @author James Till Matta
 * 
 * \tparam ParamType The floating point type returned in the uniform and non-uniform floating point distributions
 * \tparam Alpha The constant that governs the distribution, which ranges from 1/Alpha to Alpha
 * 
 * This functor class converts a floating point number drawn from the uniform
 * distribution between 0 and 1 to a floating point number draw from the
 * distribution proportional to 1/sqrt(z) on the interval 1/Alpha to Alpha
 */
template <class ParamType, double Alpha>
class GwDistribution
{
public:
    const double InverseAlpha = (1.0/Alpha);///< Constant needed for the inversion calculation
    const double QuadraticConst = (Alpha - 2.0 + InverseAlpha);///< Constant needed for the inversion calculation
    const double LinearConst = (1.0 - InverseAlpha);///< Constant needed for the inversion calculation
    /*!
     * \brief GwDistribution Constructs a GwDistribution
     */
    GwDistribution(){}
    ~GwDistribution(){}
    
    
    /*!
     * \brief operator() Returns a random sample from the distribution proportional to 1/sqrt(z)
     * \param in The floating point number drawn from the uniform distribution on [0, 1)
     * \return a floating point number draw from the distribution proportional to 1/sqrt(z) in the range 1/Alpha to Alphal
     */
    ParamType operator()(double in){return (InverseAlpha + (LinearConst*in) + (QuadraticConst*in*in));}
private:
};

}
}
#endif  //MCMC_UTILITY_MULTISAMPLER_H
