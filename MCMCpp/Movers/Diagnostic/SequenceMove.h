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
#ifndef MCMCPP_MOVERS_DIAGNOSTIC_SEQUENCEMOVER_H
#define MCMCPP_MOVERS_DIAGNOSTIC_SEQUENCEMOVER_H
// includes for C system headers
#include<stdlib.h>//needed for aligned allocation, which appears in C11 but does not appear in C++ until C++17
// includes for C++ system headers
#include<memory>//shared pointer
#include<algorithm>
#include<thread>
#include<chrono>
// includes from other libraries
// includes from MCMCpp
#include"../../Walker/Walker.h"
#include"../../Utility/ArrayDeleter.h"
#include"../../Utility/UserOjbectsTest.h"
#include"../../Utility/Misc.h"

namespace MCMC
{
namespace Mover
{

namespace Detail
{
static const int NthPrime = 50; ///<Factor to tune as a "work factor" for the calculation
}

/**
 * @class SequenceMove
 * @ingroup Movers
 * @brief An object merely increments the current walker position using a known step size for each parameter
 * @author James Till Matta
 * 
 * \tparam ParamType The floating point type to be used for the chain, float, double, long double, etc.
 * 
 * A mover that generates a monotonicly increasing sequence
 */
template <class ParamType>
class SequenceMove
{
public:
    typedef Walker::Walker<ParamType> WalkType;
    
    /*!
     * \brief SequenceMove Constructs a new sequence move object
     * \param numParams The number of parameters to work with
     * \param prngInit The seed for the random number generator
     * \param orig The original calculator class that will be copied to make the one stored internally
     */
    SequenceMove(int numParams, const ParamType* steps):
        proposal( Utility::autoAlignedAlloc<ParamType>(numParams)),
        stepSizes(Utility::autoAlignedAlloc<ParamType>(numParams), Utility::AlignedArrayDeleter<ParamType>()),
        paramCount(numParams)
    {
        for(int i=0; i<paramCount; ++i)
        {
            stepSizes.get()[i] = steps[i];
        }
    }
    
    ~SequenceMove(){Utility::delAAA(proposal);}
    
    /*!
     * \brief SequenceMove Copy constructor
     * \param rhs Original mover to be copied
     */
    SequenceMove(const SequenceMove<ParamType>& rhs):
        proposal(new ParamType[rhs.paramCount]),
        stepSizes(rhs.stepSizes),
        paramCount(rhs.paramCount){}
    
    /*!
     * \brief deleted assignment operator
     */
    SequenceMove<ParamType>& operator=(const SequenceMove<ParamType>& rhs) = delete;
    
    /*!
     * \brief setPrngSeed Sets the seed and stream number of the underlying prng
     * \param seed The seed for the prng
     * \param stream The stream number for the prng
     */
    void setPrng(long long seed, long long stream){}
    
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
        findPrimeNumber(Detail::NthPrime);
        const ParamType* currState = currWalker.getCurrState();
        //benchmarking shows that the transform was a bad choice here
        /*std::transform(currState, currState+paramCount, stepSizes.get(), proposal,
                       [](const ParamType& v1, const ParamType& v2) -> ParamType {return (v1+v2);})*/
        for(int i=0; i<paramCount; ++i)
        {
            proposal[i] = (currState[i] + stepSizes.get()[i]);
        }
        //currWalker.jumpToNewPoint(proposal, static_cast<ParamType>(0), storePoint);
        currWalker.jumpToNewPointSwap(proposal, static_cast<ParamType>(0), storePoint);
    }

    /*!
     * \brief getInitialPoints Generates an array of initial points for the walkers all set to zero
     * \param initArray The array to hold the initial points
     * \param numWalkers The number of walkers to generate initial points for
     */
    void getInitialPoints(ParamType* initArray, ParamType* auxArray, int numWalkers)
    {
        std::fill_n(initArray, numWalkers*paramCount, static_cast<ParamType>(0));
        std::fill_n(auxArray, numWalkers, static_cast<ParamType>(0));
        /*for(int j=0; j<numWalkers; ++j)
        {
            for(int i=0; i<paramCount; ++i)
            {
                initArray[j*paramCount+i] = static_cast<ParamType>(0);
            }
        }*/
    }
    
private:
    /*!
     * \brief findPrimeNumber Finds the nth prime number and returns it
     * \param n Which prime number to find
     * \return The nth prime number
     * 
     * This function exists to create a heavy CPU load, simulating a function that costs a lot to calculate
     */
    unsigned long long findPrimeNumber(int n)
    {
        int count=0;
        unsigned long long a = 2;
        while(count<n)
        {
            bool isPrime = true;
            for(unsigned long long b = 2; b * b <= a; ++b)
            {
                if(a % b == 0)
                {
                    isPrime = false;
                    break;
                }
            }
            if(isPrime) count++;
            ++a;
        }
        return (a - 1);
    }
    
    ParamType* proposal; ///Holds the new parameter set for the walker
    std::shared_ptr<ParamType> stepSizes; ///<Holds the array of step sizes, doesn't need to be replicated
    int paramCount; ///<Holds the total number of parameters
};

}

}
#endif  //MCMC_WALKER_MOVERS_DIAGNOSTIC_SEQUENCEMOVER_H
