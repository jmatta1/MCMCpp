/*!*****************************************************************************
********************************************************************************
**
** @copyright Copyright (C) 2018 James Till Matta
**
** This program is distributed in the hope that it will be useful,
** but WITHOUT ANY WARRANTY; without even the implied warranty of
** MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.
** 
********************************************************************************
*******************************************************************************/
#ifndef MCMCPP_MOVERS_DIFFERENTIALEVOLUTION_H
#define MCMCPP_MOVERS_DIFFERENTIALEVOLUTION_H
// includes for C system headers
// includes for C++ system headers
#include<memory>
// includes from other libraries
// includes from MCMCpp..
#include"../Walker/Walker.h"
#include"../Utility/MultiSampler.h"
#include"../Utility/GwDistribution.h"
#include"../Utility/UserOjbectsTest.h"
#include"../Utility/Misc.h"
#include"../Utility/ArrayDeleter.h"

namespace MCMC
{
namespace Mover
{
/**
 * @class MetropolisHastings
 * @ingroup Movers
 * @brief An object that calculates the next proposed step for a walker simply by drawing a sample from a given covariance matrix
 * @author James Till Matta
 * 
 * \tparam ParamType The floating point type to be used for the chain, float, double, long double, etc.
 * \tparam Calculator The class that calculates the log posterior and whatever else a mover may need
 * 
 * A mover that draws a random vector with mean 0 and a given covariance matrix, it adds this vector to the current
 * walker position to
 */
template <class ParamType, class Calculator>
class MetropolisHastings
{
    typedef Utility::GwDistribution<ParamType, 2, 1> DistType;
public:
    typedef Walker::Walker<ParamType> WalkType;
    static_assert(Utility::CheckCalcLogPostProb<Calculator, ParamType, ParamType*>(),
                  "MetropolisHastings: The Calculator class does not have the necessary member function with signature:\n"
                  "  'ParamType calcLogPostProb(ParamType* paramSet)'");
    static_assert(std::is_copy_constructible<Calculator>::value, "The Calculator class needs to be copy constructible.");
    
    /*!
     * @brief Constructs a new stretch move object
     * @param numParams The number of parameters to work with
     * @param prngInit The seed for the random number generator
     * @param orig The original calculator class that will be copied to make the one stored internally
     */
    MetropolisHastings(int numParams, long long prngInit, const Calculator& orig):
        decompCovMat(Utility::autoAlignedAlloc<ParamType>(numParams*numParams), Utility::AlignedArrayDeleter<ParamType>()),
        proposal(Utility::autoAlignedAlloc<ParamType>(numParams)),
        paramCount(numParams), prng(prngInit), calc(orig)
    {setIdentityCovar();}
    
    ~MetropolisHastings(){Utility::delAAA(proposal);}

    /*!
     * \brief Copy constructor
     * \param rhs Original DifferentialEvolution object to be copied
     */
    MetropolisHastings(const MetropolisHastings<ParamType, Calculator>& rhs):
        decompCovMat(rhs.decompCovMat),
        proposal(Utility::autoAlignedAlloc<ParamType>(rhs.paramCount)),
        paramCount(rhs.paramCount), diagonal(rhs.diagonal), calc(rhs.calc)
    {}
    
    /*!
     * \brief Deleted assignment operator
     */
    MetropolisHastings<ParamType, Calculator>& operator=(const MetropolisHastings<ParamType, Calculator>& rhs) = delete;
    
    
    
    /*!
     * @brief Allows the user to set a covariance matrix that is not the identity matrix
     * @param cMat the pointer to the square matrix (of size paramCount*paramCount) containing the parameter covariances
     * @return True if the matrix is accepted, False if the matrix failed certain tests or cholesky decomposition failed
     * 
     * The caller of this function retains ownership of the memory containing the covariance matrix and is responsible for deleting it
     */
    bool setCovarMat(const ParamType* cMat)
    {
        int result = testCovar(cMat);
        if(result == 1)
        {
            if(!decomposeCovarMat(cMat))
            {
                setIdentityCovar();
                return false;
            }
            else
            {
                return true;
            }
        }
        else if(result == 2)
        {
            setDiagonalDecomp(cMat);
            return true;
        }
        else
        {
            setIdentityCovar();
            return false;
        }
    }
    
    /*!
     * \brief Instructs the mover to use the simple identity covariance matrix
     */
    void setIdentityCovarMat(){setIdentityCovar();}
    
    /*!
     * \brief setPrngSeed Sets the seed and stream number of the underlying prng
     * \param seed The seed for the prng
     * \param stream The stream number for the prng
     */
    void setPrng(long long seed, long long stream){prng.setPrng(seed, stream);}
    
    /*!
     * \brief getProposal Takes the curent walker, a set of walkers to draw a target from and calculates, assumes that
     * currWalker in not in the set of walkers to select from
     * \param currWalker A reference to the current walker that we are generating a proposal for
     * \param walkerSet A pointer to the set of walkers used to generate the proposal
     * \param numWalkers The number of walkers in WalkerSet
     * \param storePoint False if the point should not be written into the chain, accepted or not
     */
    void updateWalker(WalkType& currWalker, WalkType* walkerSet, int numWalkers, bool storePoint)
    {
        //draw sample with a given covariance
        getCovarSample();
        //add the current position to the drawn sample
        const ParamType* currPos = currWalker.getCurrState();
        for(int i=0; i<paramCount; ++i)
        {
            proposal[i] += currPos[i];
        }
        ParamType currAuxVal = calc.calcLogPostProb(proposal);
        //figure out what to do with the proposal
        if((currAuxVal - currWalker.getCurrAuxData()) > (prng.getNegExponentialReal()))
        {
            currWalker.jumpToNewPointSwap(proposal, currAuxVal, storePoint);
        }
        else
        {
            currWalker.stayAtCurrentPoint(storePoint);
        }
    }

private:
    /*!
     * \brief A function to draw a random vector with the given covariance matrix, placing it in in the proposal array
     */
    void getCovarSample()
    {
        if(diagonal) getDiagonalCovarSample();
        //otherwise do the more complicated thing
        //first draw paramCount independent standard normal distributed values
        for(int i=0; i<paramCount; ++i)
        {
            proposal[i] = prng.getNormalReal();
        }
        //otherwise multiply them by the left triangular matrix from the Cholesky decomposition of the covariance matrix
        ParamType* mat = decompCovMat.get();
        for(int i=(paramCount-1); i>=0; --i)
        {//do this loop in reverse order so we do not overwrite values that we need before we are done with them
            int offset = (i*paramCount);
            /*
            //since this is a summation we *could* use kahan summation to make sure we don't lose precision for larger parameter counts
            ParamType comp = 0.0;
            ParamType sum = 0.0;
            for(int j=i; j>=0; --j)
            {
                ParamType compTerm = ((proposal[i]*mat[offset+j]) - comp);
                ParamType tempSum = (sum + compTerm);
                comp = ((tempSum - sum) - compTerm);
                sum = tempSum;
            }
            proposal[i] = sum;
            //however this is probably unnecessary for the approximations inherent in an MCMC scheme anyways
            */
            proposal[i] *= mat[offset+i];
            for(int j=(i-1); j>=0; --j)
            {
                proposal[i] += (proposal[j]*mat[offset+j]);
            }
        }
    }

    /*!
     * @brief Calculates the sample in the simple case of a purely diagonal sample matrix, storing it in the proposal array
     */
    void getDiagonalCovarSample()
    {
        //draw the normals multiplied by the diagonal element and we are done
        ParamType* mat = decompCovMat.get();
        for(int i=0; i<paramCount; ++i)
        {
            proposal[i] = (prng.getNormalReal()*mat[paramCount*i+i]);
        }
    }
    
    /*!
     * \brief Tests the covariance matrix for suitability
     * \param cMat The covariance matrix to be tested.
     * \return 0 for a matrix that is not symmetrix or has a non-positive value on the diagonal, 2 for a diagonal matrix that has only positive values, 1 in all other cases
     */
    int testCovar(const ParamType* cMat)
    {
        bool isDiagonal = true;
        for(int i=1; i<paramCount; ++i)
        {
            //test that all the diagonals are greater than zero
            if(cMat[i*paramCount+i] <= 0.0) return 0;//zero or negative value in diagonal
            //test that it is symmetric
            for(int j=0; j<i; ++j)
            {
                ParamType val = cMat[i*paramCount+j];
                if(val != cMat[j*paramCount+i]) return 0;//not symmetric
                if(val != 0.0) isDiagonal = false;
            }
        }
        //check if we discovered that it is diagonal
        if(isDiagonal) return 2;
        //it is not diagonal but passes the other tests
        return 1;
    }
    
    /*!
     * @brief Perform a Cholesky decomposition of the (tested) covariance matrix and store the result
     * @param cMat The covariance matrix to be decomposed
     * @return True for successful decomposition, false for failure
     */
    bool decomposeCovarMat(const ParamType* cMat)
    {
        ParamType* temp = decompCovMat.get();
        for(int i=0; i<paramCount; ++i)
        {
            for(int j=0; j<paramCount; ++j)
            {
                if(i==j)//diagonal term case
                {
                    ParamType value = cMat[i*paramCount+j];
                    for(int k=0; k<j; ++k)
                    {
                        value -= (temp[j*paramCount+k]*temp[j*paramCount+k]);
                    }
                    if(value >= 0.0)
                    {
                        temp[i*paramCount+j] = std::sqrt(value);
                    }
                    else//something broke return that the matrix was not valid
                    {
                        return false;
                    }
                }
                else if(j<i)//non diagonal case with non-zero values
                {
                    ParamType value = cMat[i*paramCount+j];
                    for(int k=0; k<j; ++k)
                    {
                        value -= (temp[i*paramCount+k]*temp[j*paramCount+k]);
                    }
                    temp[i*paramCount+j] = (value/temp[j*paramCount+j]);
                }
                else//non diagonal case were the terms are zero
                {
                    temp[i*paramCount+j] = 0.0;
                }
            }
        }
        diagonal = false;
        return true;
    }
    
    /*!
     * \brief When the covaraiance matrix is diagonal, this merely sets the decomposed matrix to the square roots of the covar matrices elements
     * \param cMat the diagonal covariance matrix
     */
    void setDiagonalDecomp(const ParamType* cMat)
    {
        ParamType* temp = decompCovMat.get();
        for(int i=0; i<paramCount; ++i)
        {
            int offset = i*paramCount;
            for(int j=0; j<paramCount; ++j)
            {
                if(i==j)
                {
                    temp[offset + j] = std::sqrt(cMat[offset + j]);
                }
                else
                {
                    temp[offset + j] = 0.0;
                }
            }
        }
        diagonal = true;
    }
    
    /*!
     * \brief Sets the decomposition matrix to the identity matrix (used as a default for bad covar matrices)
     */
    void setIdentityCovar()
    {
        ParamType* temp = decompCovMat.get();
        for(int i=0; i<paramCount; ++i)
        {
            int offset = i*paramCount;
            for(int j=0; j<paramCount; ++j)
            {
                if(i==j)
                {
                    temp[offset + j] = 1.0;
                }
                else
                {
                    temp[offset + j] = 0.0;
                }
            }
        }
        diagonal = true;
    }
    
    std::shared_ptr<ParamType> decompCovMat;
    ParamType* proposal = nullptr;
    int paramCount;
    bool diagonal = true;
    Utility::MultiSampler<ParamType, DistType> prng;
    Calculator calc;
};
}
}
#endif  //MCMCPP_MOVERS_DIFFERENTIALEVOLUTION_H
