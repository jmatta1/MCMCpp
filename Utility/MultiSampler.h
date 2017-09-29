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
 * \tparam CustomDistribution The class that performs the cdf inversion to generate a custom distributed floating point
 * \tparam ParamType The floating point type returned in the uniform and non-uniform floating point distributions
 * 
 * This class manages the c++ PRNG underlying these and draws samples as needed from that PRNG
 */
template <class CustomDistribution, class ParamType>
class MultiSampler
{
public:
    /*!
     * \brief MultiSampler Constructs a multisampler with a given seed
     * \param seed The 64 bit integer that will serve the engine as a seed
     */
    MultiSampler(long long seed):engine(seed),realDist(){}
    ~MultiSampler(){}
    
    /*!
     * \brief getCustomSample returns a sample from the given probability distribution
     * \return ParamType floating point number drawn from the given distribution
     */
    ParamType getCustomSample(){return customDist(realDist(engine));}
    /*!
     * \brief getUniformReal returns a sample from the flat probability distribution between 0 and 1
     * \return ParamType floating point draw from the flat distribution between 0 and 1
     */
    ParamType getUniformReal(){return realDist(engine);}
    
    /*!
     * \brief getNonOffSetInt Gets a random integer in the range [0, max)
     * \param max The maximum value of the range plus 1
     * \return A pseudo random integer in the range [0, max)
     */
    int getNonOffSetInt(int max){return (engine()%max);}
    /*!
     * \brief getNonOffSetInt Gets a random integer in the range [offset, max)
     * \param max The maximum value of the range plus 1
     * \return A pseudo random integer in the range [offset, max)
     */
    int getOffSetInt(int offset, int max){return ((engine()%(max-offset))+offset);}
    /*!
     * \brief getNonOffSetInt Gets a random 64-bit integer in the range [0, max)
     * \param max The maximum value of the range plus 1
     * \return A pseudo random 64-bit integer in the range [0, max)
     */
    long long getNonOffSetLongLong(long long max){return (engine()%max);}
    /*!
     * \brief getNonOffSetInt Gets a random 64-bit integer in the range [offset, max)
     * \param max The maximum value of the range plus 1
     * \return A pseudo random 64-bit integer in the range [offset, max)
     */
    long long getOffSetLongLong(long long offset, long long max){return ((engine()%(max-offset))+offset);}
private:
    std::mt19937_64 engine;
    std::uniform_real_distribution<ParamType> realDist;
    CustomDistribution customDist;
};

}
}
#endif  //MCMC_UTILITY_MULTISAMPLER_H
