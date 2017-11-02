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
#include<random>
// includes from other libraries
// includes from MCMC

namespace MarkovChainMonteCarlo
{

namespace Utility
{

/**
 * @class MultiSampler
 * @ingroup Utility
 * @brief The class used to generate all kinds of pseudo random numbers needed for MCMC, ints, uniform floating point [0...1], non-uniform floating point
 * @author James Till Matta
 * 
 * \tparam ParamType The floating point type returned in the uniform and non-uniform floating point distributions
 * \tparam CustomDistribution The class that performs the cdf inversion to generate a custom distributed floating point
 * 
 * This class manages the c++ PRNG underlying these and draws samples as needed from that PRNG
 */
template <class ParamType, class CustomDistribution>
class MultiSampler
{
public:
    /*!
     * \brief MultiSampler Constructs a multisampler with a given seed
     * \param seed The 64 bit integer that will serve the engine as a seed
     */
    MultiSampler(long long seed=0):engine(seed),realDist(), normDist(){}
    ~MultiSampler(){}
    
    /*!
     * \brief setPrngSeed Sets the random number engines seed to see
     * \param seed The new seed for the random number engine
     */
    inline void setPrngSeed(long long seed){engine.seed(seed);}
    
    /*!
     * \brief getCustomSample returns a sample from the given probability distribution
     * \return ParamType floating point number drawn from the given distribution
     */
    inline ParamType getCustomSample(){return customDist(realDist(engine));}
    /*!
     * \brief getUniformReal returns a sample from the flat probability distribution between 0 and 1
     * \return ParamType floating point drawn from the flat distribution between 0 and 1
     */
    inline ParamType getUniformReal(){return realDist(engine);}
    
    /*!
     * \brief getExponentialReal returns a sample from the exponential probability distribution with parameter 1
     * \return ParamType floating point drawn from the exponential distribution
     */
    inline ParamType getExponentialReal(){return expDist(engine);}
    
    /*!
     * \brief getNegExponentialReal returns the negative of asample from the exponential probability distribution with parameter 1
     * \return ParamType floating point drawn from the exponential distribution, then made negative, which gives it the same shape as Log(x) where x is uniform on [0,1)
     */
    inline ParamType getNegExponentialReal(){return -expDist(engine);}
    
    /*!
     * \brief getNormalReal returns a sample from the normal probability distribution centered at 0 and with variance 1
     * \return ParamType floating point drawn from the normal distribution centered at 0 and with variance 1
     */
    inline ParamType getNormalReal(){return normDist(engine);}
    
    /*!
     * \brief getNonOffSetInt Gets a random integer in the range [0, max)
     * \param max The maximum value of the range plus 1
     * \return A pseudo random integer in the range [0, max)
     */
    inline int getNonOffSetInt(int max){return (engine()%max);}
    /*!
     * \brief getNonOffSetInt Gets a random integer in the range [offset, max)
     * \param max The maximum value of the range plus 1
     * \return A pseudo random integer in the range [offset, max)
     */
    inline int getOffSetInt(int offset, int max){return ((engine()%(max-offset))+offset);}
    /*!
     * \brief getNonOffSetInt Gets a random 64-bit integer in the range [0, max)
     * \param max The maximum value of the range plus 1
     * \return A pseudo random 64-bit integer in the range [0, max)
     */
    inline long long getNonOffSetLongLong(long long max){return (engine()%max);}
    /*!
     * \brief getNonOffSetInt Gets a random 64-bit integer in the range [offset, max)
     * \param max The maximum value of the range plus 1
     * \return A pseudo random 64-bit integer in the range [offset, max)
     */
    inline long long getOffSetLongLong(long long offset, long long max){return ((engine()%(max-offset))+offset);}
private:
    std::mt19937_64 engine;///<The base random number generator engine that is used for everything
    std::uniform_real_distribution<ParamType> realDist;///<The adapter that gives uniform real numbers between 0 and 1
    std::exponential_distribution<ParamType> expDist;///<This adapter gives an exponential distribution with form e^-x, which is the shape of -log(x) when x is a uniform distribution
    std::normal_distribution<ParamType> normDist;///<The adapter that gives normally distributed real numbers with mean 0 and variance 1
    CustomDistribution customDist;///<The CDF inversion class that takes a number from 0 to 1 and converts it to a custom distribution
};

}
}
#endif  //MCMC_UTILITY_MULTISAMPLER_H
