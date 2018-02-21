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
#ifndef MCMCPP_ANALYSIS_COVARIANCEMATRIX_H
#define MCMCPP_ANALYSIS_COVARIANCEMATRIX_H
// includes for C system headers
// includes for C++ system headers
#include<cmath>//for sqrt
#include<cassert>//for assert
// includes from other libraries
// includes from MCMCpp
#include"../Utility/Misc.h"
#include"../Chain/ChainStepIterator.h" //the iterator to step through chains with


namespace MCMC
{
namespace Analysis
{

/*!
 * @class CovarianceMatrix
 * @ingroup Analysis
 * @brief A class to calculate the covariance matrix of a given set of MCMC samples
 * @author James Till Matta
 * 
 * @tparam ParamType The floating point type for which this calculation will be carried out
 * 
 * This class offers two versions of the covariance calculation. In the first,
 * every sample is used and possible effects of autocorrelation are ignored. In
 * the second only a subset of samples are used according to some given slicing
 * (only use every nth sample) a crude method for accounting for
 * autocorrelation.
 */
template<class ParamType>
class CovarianceMatrix
{
public:
    typedef Chain::ChainStepIterator<ParamType> IttType;
    
    /*!
     * \brief Constructor for the Covariance Matrix Calculator
     * \param numParams The number of parameters in the matrix
     * \param numWalkers The number of walkers that are in ensemble
     */
    CovarianceMatrix(int numParams, int numWalkers);
    
    /*!
     * \brief Destructor which uses the deleteAAA function to free all allocated arrays
     */
    ~CovarianceMatrix();
    
    /*!
     * \brief Uses every sliceInterval'th sample to calculate the covariance of the markov chains
     * \param start The start chain step iterator
     * \param end The end chain step iterator
     * \param sliceInterval The relative index of sample to use, must be >= 1. a value of 1 uses every sample, 2 uses every other sample, so on and so forth
     */
    void calculateCovar(IttType start, IttType end, int sliceInterval=1);
    
    /*!
     * \brief Gets the specified element of the covariance matrix
     * \param row The row of the element to select
     * \param col The column of the element to select
     * \return The desired element of the covariance matrix
     */
    ParamType getCovarianceMatrixElement(int row, int col){return covarMat[row*pCount+col];}
    
    /*!
     * \brief Gets the specified element of the correlation matrix
     * \param row The row of the element to select
     * \param col The column of the element to select
     * \return The desired element of the correlation matrix
     */
    ParamType getCorrelationMatrixElement(int row, int col){return corrMat[row*pCount+col];}
    
private:
    /*!
     * \brief Takes the completed calculation of pProducts and makes the averages and covariance matrix
     * \param numUsedSteps The number of steps that were used in the calculation of the averages and the parameter products
     */
    void finalizeMatrix(int numUsedSteps);
    
    /*!
     * \brief Takes a point from a markov chain to process into the covariance matrix
     * \param pt A pointer to the start of the point in question
     */
    void processPoint(ParamType* pt);
    
    /*!
     * \brief Adds a product to the appropriate slots of the covariance matrix
     * \param i Index of the first parameter in the product
     * \param j Index of the second parameter in the product
     * \param val The product of the two parameters
     */
    void addProduct(int i, int j, ParamType val);
    
    /*!
     * \brief Accumulates val in the i'th index of pAvgs array using kahan summation
     * \param i The index of the parameter
     * \param val The value of the parameter in the chain
     */
    void addAverage(int i, ParamType val);
    
    /*!
     * \brief Zeros the arrays and their compensation values
     */
    void zeroStorage();
    
    ParamType* covarMat = nullptr; ///<Holds the final version of the covariance matrix
    ParamType* corrMat = nullptr; ///<Holds the correlation matrix (a normalized version of the covariance matrix)
    ParamType* pProducts = nullptr; ///<Holds the Xi*Yi part of the covariance matrix formula
    ParamType* pProductsComp = nullptr; ///<Holds the kahan summation compensations for the pProducts
    ParamType* pAvgs = nullptr; ///<Holds the sums for parameter averaging (and then the average itself)
    ParamType* pAvgsComp = nullptr; ///<Holds the kahan summation compensations for the parameter averages
    int pCount = 0; ///<Holds the parameter count
    int wCount = 0; ///<Holds the walker count
};

template<class ParamType>
CovarianceMatrix<ParamType>::CovarianceMatrix(int numParams, int numWalkers):
    covarMat(Utility::autoAlignedAlloc<ParamType>(numParams*numParams)),
    corrMat(Utility::autoAlignedAlloc<ParamType>(numParams*numParams)),
    pProducts(Utility::autoAlignedAlloc<ParamType>(numParams*numParams)),
    pProductsComp(Utility::autoAlignedAlloc<ParamType>(numParams*numParams)),
    pAvgs(Utility::autoAlignedAlloc<ParamType>(numParams*numParams)),
    pAvgsComp(Utility::autoAlignedAlloc<ParamType>(numParams*numParams)),
    pCount(numParams), wCount(numWalkers)
{
    assert(pCount > 0);
    assert(wCount > 0);
}

template<class ParamType>
CovarianceMatrix<ParamType>::~CovarianceMatrix()
{
    Utility::delAAA(covarMat);
    Utility::delAAA(corrMat);
    Utility::delAAA(pProducts);
    Utility::delAAA(pProductsComp);
    Utility::delAAA(pAvgs);
    Utility::delAAA(pAvgsComp);
}

template<class ParamType>
void CovarianceMatrix<ParamType>::calculateCovar(IttType start, IttType end, int sliceInterval=1)
{
    zeroStorage();
    int index = 0;
    int numPts = 0;
    for(IttType itt(start); itt!=end; ++itt)
    {
        if((index % sliceInterval) != 0)
        {
            ++index;
            continue; //skip samples that are not the sliceInterval'th sample
        }
        for(int i=0; i<wCount; ++i)
        {
            processPoint( ((*itt)+(i*pCount)) );
        }
        ++index;
        ++numPts;
    }
    finalizeMatrix(numPts);
}

template<class ParamType>
void CovarianceMatrix<ParamType>::calculateCovarSlicing(IttType start, IttType end, int sliceInterval=1)
{
    zeroStorage();
    int index = 0;
    int numPts = 0;
    for(IttType itt(start); itt!=end; ++itt)
    {
        if((index % sliceInterval) != 0)
        {
            ++index;
            continue; //skip samples that are not the sliceInterval'th sample
        }
        for(int i=0; i<wCount; ++i)
        {
            processPoint( ((*itt)+(i*pCount)) );
        }
        ++index;
        ++numPts;
    }
    finalizeMatrix(numPts);
}

template<class ParamType>
void CovarianceMatrix<ParamType>::finalizeMatrix(int numUsedSteps)
{
    //first get the actual n used
    ParamType ptCount = static_cast<ParamType>(numUsedSteps*wCount);
    //then divide out the averages
    for(int i=0; i<pCount; ++i)
    {
        pAvgs[i] /= ptCount;
    }
    //then calculate the elements of the covaraiance matrix
    for(int i=0; i<pCount; ++i)
    {
        for(int j=0; j<=i; ++j)
        {
            ParamType temp = ((pProducts[i*pCount+j]/ptCount) - (pAvgs[i]*pAvgs[j]));
            covarMat[i*pCount+j] = temp;
            covarMat[j*pCount+i] = temp;
        }
    }
    //then calculate the elements of the correlation matrix
    //first iterate row by row, get the sqrt(variance) of the diagonal of that row and divide everything in the row by it
    for(int i=0; i<pCount; ++i)
    {
        ParamType sqrtVariance = std::sqrt(covarMat[i*pCount+i]);
        for(int j=0; j<pCount; ++j)
        {
            corrMat[i*pCount+j] = (covarMat[i*pCount+j]/sqrtVariance);
        }
    }
    //then iterate column by colum, get the sqrt(variance) of the diagonal of that column and divide everyint in the column by it
    for(int j=0; j<pCount; ++j)
    {
        ParamType sqrtVariance = std::sqrt(covarMat[j*pCount+j]);
        for(int i=0; i<pCount; ++i)
        {
            corrMat[i*pCount+j] /= sqrtVariance;
        }
    }
}

template<class ParamType>
void CovarianceMatrix<ParamType>::processPoint(ParamType* pt)
{
    //cycle through value pairs
    for(int i=0; i<pCount; ++i)
    {
        //before handling the full pair, add the value to the averages
        addAverage(i, pt[i]);
        for(int j=0; j<=i; ++j)
        {//handle the full pair
            addProduct(i,j, (pt[i]*pt[j]));
        }
    }
}

template<class ParamType>
void CovarianceMatrix<ParamType>::addProduct(int i, int j, ParamType val)
{
    //first handle the 'natural' order
    int offset = i*pCount;
    ParamType y = (val - pProductsComp[offset+j]);
    ParamType t = (pProducts[offset+j] + y);
    pProductsComp[offset+j] = ((t-pProducts[offset+j])-y);
    pProducts[offset+j] = t;
    
    //then, if not on a diagonal handle the opposite order
    if(i!=j)
    {
        offset = j*pCount;
        y = (val - pProductsComp[offset+i]);
        t = (pProducts[offset+i] + y);
        pProductsComp[offset+i] = ((t-pProducts[offset+i])-y);
        pProducts[offset+i] = t;
    }
}

template<class ParamType>
void CovarianceMatrix<ParamType>::addAverage(int i, ParamType val)
{
    ParamType y = (val - pAvgsComp[i]);
    ParamType t = (pAvgs[i] + y);
    pAvgsComp[i] = ((t-pAvgs[i])-y);
    pAvgs[i] = t;
}

template<class ParamType>
void CovarianceMatrix<ParamType>::zeroStorage()
{
    ParamType zero = static_cast<ParamType>(0);
    for(int i=0; i<pCount; ++i)
    {
        pAvgs[i] = zero;
        pAvgsComp[i] = zero;
        int offset = i*pCount;
        for(int j=0; j<pCount; ++j)
        {
            covarMat[offset + j] = zero;
            corrMat[offset + j] = zero;
            pProducts[offset + j] = zero;
            pProductsComp[offset + j] = zero;
        }
    }
}

}
}
#endif  //MCMCPP_ANALYSIS_COVARIANCEMATRIX_H
