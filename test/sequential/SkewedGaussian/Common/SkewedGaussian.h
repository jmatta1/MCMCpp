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
#ifndef COMMON_SKEWEDGAUSSIAN_H
#define COMMON_SKEWEDGAUSSIAN_H
// includes for C system headers
// includes for C++ system headers
// includes from other libraries
// includes from MCMC

/**
 * @class SkewedGaussianTwoDim
 * @brief A parameterized skewed gaussian in two dimensions
 * @author James Till Matta
 * 
 * \tparam ParamType The floating point type to be used for the chain, float, double, long double, etc.
 */
template <class ParamType>
class SkewedGaussianTwoDim
{
public:
    /*!
     * \brief SkewedGaussianTwoDim Constructor for the Two Dimensional Skewed Gaussian likelihood function
     * \param eps The skew parameter for the distribution
     */
    SkewedGaussianTwoDim(const ParamType& eps):epsilon(eps){}
    
    /*!
     * \brief ~SkewedGaussianTwoDim simple destructor
     */
    ~SkewedGaussianTwoDim(){}
    
    /*!
     * \brief SkewedGaussianTwoDim Copy Constructor for the Two Dimensional Skewed Gaussian likelihood function
     * \param rhs The original object to be copied
     */
    SkewedGaussianTwoDim(const SkewedGaussianTwoDim& rhs):epsilon(rhs.epsilon){}
    
    /*!
     * \brief calcLogPostProb Calculates the log of the posterior probability distribution
     * \param paramSet A pointer to the set of parameters
     * \return The log of the posterior probability distribution
     */
    ParamType calcLogPostProb(ParamType* paramSet)
    {
        ParamType temp1 = (paramSet[0] - paramSet[1]);
        ParamType temp2 = (paramSet[0] + paramSet[1]);
        return (((temp1*temp1)/epsilon + (temp2*temp2))/static_cast<ParamType>(-2));
    }

private:
    ParamType epsilon; ///<Stores the skew parameter of the distribution
};


#endif  //COMMON_SKEWEDGAUSSIAN_H
