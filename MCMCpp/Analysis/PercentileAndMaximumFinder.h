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
#ifndef MCMCPP_ANALYSIS_PERCENTILES_H
#define MCMCPP_ANALYSIS_PERCENTILES_H
// includes for C system headers
// includes for C++ system headers
#include<limits>//for numeric_limits
// includes from other libraries
// includes from MCMCpp
#include"../Utility/Misc.h"
#include"../Chain/ChainStepIterator.h" //the iterator to step through chains with

namespace MCMC
{
namespace Analysis
{

/*!
 * @class PercentileAndMaximumFinder
 * @ingroup Analysis
 * @brief A class to calculate the single and bi-axial histograms of the points in chains
 * @author James Till Matta
 * 
 * @tparam ParamType The floating point type for which this calculation will be carried out
 * 
 * A method to do this could work by storing all the sets of parameters in sorted lists.
 * This method would take size proportional to the size of the chain, which could be very
 * large. Instead this class uses the somewhat less accurate method of binning the parameter
 * values and finding percentiles from the cumulative sums of those histograms.
 * The user can skip parts of the chain to use every nth step in or simply use every point.
 * Choosing a reasonable binning is not as memory critical here as it is for CornerHistograms.
 * The size taken will only grow as 4*n*b where n is the number of parameters and b is
 * the number of bins. However, for finding the distribution maximum a binning that is too
 * fine could make the determined maximum unreliable (because the histogram could
 * become jagged). The maximum finding algorithm will be limitted in accuracy by
 * the binning applied and the range of values of a given parameter
 */
template<class ParamType>
class PercentileAndMaximumFinder
{
public:
    typedef Chain::ChainStepIterator<ParamType> IttType;
    /*!
     * \brief Constructor for Percentile and maximum finder
     * \param numParams The number of parameters per point
     * \param numWalkers The number of walkers per step in the chain
     * \param binsPerAxis The desired number of bins per axis in the histograms, 1000 bins per axis is a somewhat sane default value since we have no 2D histograms
     */
    PercentileAndMaximumFinder(int numParams, int numWalkers, int binsPerAxis=1000);

    ~PercentileAndMaximumFinder(){delete[] hists; delete[] cumSum; delete[] paramBounds;}
    
    /*!
     * \brief Histograms the points in the chain into 1D histograms to use for finding percentiles
     * \param start The step iterator at the start of the chain
     * \param end The step iterator at the end of the chain
     * \param sliceInterval The relative index of sample to use, must be >= 1. a value of 1 uses every sample, 2 uses every other sample, so on and so forth
     */
    void calculateHistograms(IttType start, IttType end, int sliceInterval=1);
    
    /*!
     * \brief Takes a given parameter and value and finds what the percentile for that value is
     * \param pIndex The index of the parameter to find the percentile for
     * \param val The value for which a percentile desired
     * \return The percentile of the value passed if it is within range, otherwise, -1
     */
    ParamType getPercentileFromValue(int pIndex, const ParamType& val);
    
    /*!
     * \brief Takes a given percentile and finds the value corresponding to that percentile
     * \param pIndex The index of the parameter to find the value for
     * \param val The percentile desired
     * \return The value at the percentile passed if percentile is on [0,100], otherise the minumum value of the range is returned
     */
    ParamType getValueFromPercentile(int pIndex, const ParamType& per);
    
    /*!
     * \brief Finds the center of the most probable bin for a given parameter
     * \param pIndex The parameter to find the maximum for
     * \return The center of the most full bin
     */
    ParamType getValueOfPeak(int pIndex);
    
    /*!
     * \brief Gives the minimum value of the parameter's range
     * \param pIndex The parameter whose minimum is to be given
     * \return The bottom of a parameter's range of values
     */
    ParamType getParamMinimum(int pIndex){return paramBounds[2*pIndex];}
    
    /*!
     * \brief Gives the maximum value of the parameter's range
     * \param pIndex The parameter whose maximum is to be given
     * \return The top of a parameter's range of values
     */
    ParamType getParamMaximum(int pIndex){return (paramBounds[2*pIndex]+(paramBounds[2*pIndex+1]*bCount));}

private:
    /*!
     * \brief Calculates the cumulative sums from the histograms
     */
    void calcCumSums();
    
    /*!
     * \brief Takes a point and increments the histograms as needed using the point
     * \param pt The point to be used to increment histograms
     */
    void handlePoint(ParamType* pt);
    
    /*!
     * \brief Finds the minimum and maximum value of every parameter
     * \param start The iterator at the start of the chain
     * \param end The iterator at the end of the chain
     * \param sliceInterval The relative index of sample to use, must be >= 1. a value of 1 uses every sample, 2 uses every other sample, so on and so forth
     */
    void findBinning(IttType start, IttType end, int sliceInterval);
    
    /*!
     * \brief Checks if each parameter constitutes an extreme (min or max) for a given point
     * \param pt The point whose values are to be tested
     */
    void testPtExtremes(ParamType* pt);
    
    /*!
     * \brief Resets all the storage in the class to zero, or relevant value for starting
     */
    void zeroStorage();
    
    /*!
     * \brief Returns the sign of val, -1 for negative numbers 1 for positive, 0 for zeros
     * \param val The value whose sign is to be determined
     * \return -1 if val is negative, 1 if val is positive, 0 if val is 0
     */
    int sign(ParamType val)
    {
        return (ParamType(0) < val) - (val < ParamType(0));
    }
    
    const ParamType expandFraction = static_cast<ParamType>(1.001); ///<Holds the margin for expanding negative bounds downward and positive bounds upward
    const ParamType contractFraction = static_cast<ParamType>(0.999);///<Holds the margin for expanding negative bounds *upward* and positive bounds downward
    const ParamType minSize = static_cast<ParamType>(0.001);///<Holds the minimum magnitude of a range
    
    int* hists = nullptr;///<Holds a histogram of each parameter
    int* cumSum = nullptr;///<Holds a cumulative sum of each histogram
    ParamType* paramBounds = nullptr; ///<Holds the min and max of every parameter
    int numPoints;///<Holds the total number of points that went into generating the histograms
    int pCount;///<Holds the number of parameters per point
    int wCount;///<Holds the number of walkers per step
    int bCount;///<Holds the number of bins per histogram
    int cbCount;///<Holds the number of bins per cumulative sum histogram
    bool binned = false;///<Holds whether or not the user has supplied data to generate the histograms from
};

template<class ParamType>
PercentileAndMaximumFinder<ParamType>::PercentileAndMaximumFinder(int numParams, int numWalkers, int binsPerAxis):
    pCount(numParams), wCount(numWalkers), bCount(binsPerAxis), cbCount(binsPerAxis+1)
{
    //allocate the array of hists
    hists = new int[bCount*pCount];
    //allocate the array of cumulative sums
    cumSum = new int[cbCount*pCount];
    //allocate the set of parameter bounds and set them to extreme values
    paramBounds = new ParamType[2*pCount];
}

template<class ParamType>
void PercentileAndMaximumFinder<ParamType>::calculateHistograms(IttType start, IttType end, int sliceInterval)
{
    //reset our storage
    zeroStorage();
    //find the binning parameters
    findBinning(start, end, sliceInterval);
    //iterate through the chain steps
    int index = 0;
    for(IttType itt(start); itt!=end; ++itt)
    {
        if((index % sliceInterval) != 0)
        {
            ++index;
            continue; //skip samples that are not the sliceInterval'th sample
        }
        //now iterate through every point
        for(int i=0; i<wCount; ++i)
        {
            handlePoint( ((*itt)+(i*pCount)) );
        }
        ++index;
    }
    //calculate the cumulative sums
    calcCumSums();
    binned = true;
}

template<class ParamType>
ParamType PercentileAndMaximumFinder<ParamType>::getPercentileFromValue(int pIndex, const ParamType& val)
{
    
}

template<class ParamType>
ParamType PercentileAndMaximumFinder<ParamType>::getValueFromPercentile(int pIndex, const ParamType& per)
{
    
}

template<class ParamType>
ParamType PercentileAndMaximumFinder<ParamType>::getValueOfPeak(int pIndex)
{
    
}

template<class ParamType>
void PercentileAndMaximumFinder<ParamType>::calcCumSums()
{
    for(int i=0; i<pCount; ++i)
    {
        for(int j=0; j<bCount; ++j)
        {
            cumSum[i*cbCount+j+1] = (cumSum[i*cbCount+j] + hists[i*bCount+j]);
        }
    }
}

template<class ParamType>
void PercentileAndMaximumFinder<ParamType>::handlePoint(ParamType* pt)
{
    //iterate through the parameters
    for(int i=0; i<pCount; ++i)
    {
        //first handle the single parameter histogram
        int iBin = ((pt[i] - paramBounds[2*i])/paramBounds[2*i+1]);
        ++(hists[i*bCount+iBin]);
    }
}

template<class ParamType>
void PercentileAndMaximumFinder<ParamType>::findBinning(IttType start, IttType end, int sliceInterval)
{
    //first find the extrema of the points to use
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
            testPtExtremes( ((*itt)+(i*pCount)) );
        }
        ++index;
        ++numPts;
    }
    //now tweak the parameter bounds slightly, preventing identical values etc
    for(int i=0; i<pCount; ++i)
    {
        //check for degeneracy
        if(paramBounds[2*i] == paramBounds[2*i+1])
        {
            if(paramBounds[2*i] != static_cast<ParamType>(0))
            {
                if(sign(paramBounds[2*i]) == 1)//positive value
                {paramBounds[2*i] *= contractFraction; paramBounds[2*i+1] *= expandFraction;}
                else
                {paramBounds[2*i] *= expandFraction; paramBounds[2*i+1] *= contractFraction;}
            }
            else
            {
                paramBounds[2*i] = -minSize;
                paramBounds[2*i+1] = minSize;
            }
        }
        else
        {//now tweak the bounds slightly so they are inclusive
            int ans = sign(paramBounds[2*i]);
            if(ans == -1){paramBounds[2*i] *= expandFraction;}
            else if(ans == 0){paramBounds[2*i] = -minSize;}
            else{paramBounds[2*i] *= contractFraction;}
            ans = sign(paramBounds[2*i+1]);
            if(ans == -1){paramBounds[2*i+1] *= expandFraction;}
            else if(ans == 0){paramBounds[2*i+1] = minSize;}
            else{paramBounds[2*i+1] *= contractFraction;}
        }
    }
    //now use those parameter bounds to calculate the width per bin storing it in the no longer needed upper bound
    for(int i=0; i<pCount; ++i)
    {
        paramBounds[2*i+1] = ((paramBounds[2*i+1]-paramBounds[2*i])/static_cast<ParamType>(bCount));
    }
}

template<class ParamType>
void PercentileAndMaximumFinder<ParamType>::testPtExtremes(ParamType* pt)
{
    for(int i=0; i<pCount; ++i)
    {
        if(pt[i] < paramBounds[2*i])
        {
            paramBounds[2*i] = pt[i];
        }
        if(pt[i] > paramBounds[2*i+1])
        {
            paramBounds[2*i+1] = pt[i];
        }
    }
}

template<class ParamType>
void PercentileAndMaximumFinder<ParamType>::zeroStorage()
{
    //zero all the histograms
    int arrSize = bCount*pCount;
    for(int i=0; i<arrSize; ++i)
    {
        hists[i] = 0;
    }
    //zero the cumulative sum
    arrSize = cbCount*pCount;
    for(int i=0; i<arrSize; ++i)
    {
        cumSum[i] = 0;
    }
    //set the parameter bounds to their initial positions
    for(int i=0; i<pCount; ++i)
    {
        paramBounds[2*i] = std::numeric_limits<ParamType>::max();
        paramBounds[2*i+1] = std::numeric_limits<ParamType>::min();
    }
}

}
}
