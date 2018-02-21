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
#include<cassert>//for assert
#include<sstream>//for generating file names for saving histograms
#include<fstream>//for saving histograms
#include<string>//for passing strings around for file names
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

    /*!
     * \brief Destructor
     */
    ~PercentileAndMaximumFinder();
    
    /*!
     * \brief Histograms the points in the chain into 1D histograms to use for finding percentiles
     * \param start The step iterator at the start of the chain
     * \param end The step iterator at the end of the chain
     * \param sliceInterval The relative index of sample to use, must be >= 1. a value of 1 uses every sample, 2 uses every other sample, so on and so forth
     */
    void processChainData(IttType start, IttType end, int sliceInterval=1);
    
    /*!
     * \brief Takes a given parameter and value and finds what the percentile for that value is
     * \param pIndex The index of the parameter to find the percentile for
     * \param val The value for which a percentile desired
     * \return The percentile of the value passed (in the range [0,100]) if it is within range, otherwise, -1
     */
    ParamType getPercentileFromValue(int pIndex, const ParamType& val);
    
    /*!
     * \brief Takes a given percentile and finds the value corresponding to that percentile
     * \param pIndex The index of the parameter to find the value for
     * \param val The percentile desired (in the range [0,100])
     * \return The value at the percentile passed if percentile is on [0,100], otherise the minumum value of the range minus 1e4 is returned
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
    
    /*!
     * \brief Writes the histograms, normal and cumulative sums, to files whose names start with fileNameBase
     * \param fileNameBase The main component of the file name to which sub identifiers are prefixed on
     */
    void writeHistogramsInCsvFormat(const std::string& fileNameBase);

private:
    /*!
     * \brief Writes the specified histogram into csv format
     * \param histNum The parameter number of the histogram
     * \param fileNameBase The base filename from which the output filename will be calculated
     */
    void write1dHistCsv(int histNum, const std::string& fileNameBase);
    
    /*!
     * \brief Writes the specified cumulative sum into csv format
     * \param histNum The parameter number of the cumulative sum
     * \param fileNameBase The base filename from which the output filename will be calculated
     */
    void write1dCumSumCsv(int histNum, const std::string& fileNameBase);
    
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
    int sign(ParamType val);
    
    const ParamType expandFraction = static_cast<ParamType>(1.001); ///<Holds the margin for expanding negative bounds downward and positive bounds upward
    const ParamType contractFraction = static_cast<ParamType>(0.999);///<Holds the margin for expanding negative bounds *upward* and positive bounds downward
    const ParamType minSize = static_cast<ParamType>(0.001);///<Holds the minimum magnitude of a range
    
    int* hists = nullptr;///<Holds a histogram of each parameter
    int* cumSum = nullptr;///<Holds a cumulative sum of each histogram
    ParamType* paramBounds = nullptr; ///<Holds the min and max of every parameter
    int numPoints=0;///<Holds the total number of points that went into generating the histograms
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
    assert(pCount > 0);
    assert(wCount > 0);
    assert(bCount > 1);
    //allocate the array of hists
    hists = new int[bCount*pCount];
    //allocate the array of cumulative sums
    cumSum = new int[cbCount*pCount];
    //allocate the set of parameter bounds and set them to extreme values
    paramBounds = new ParamType[2*pCount];
    //preemptively zero everything
    zeroStorage();
}

template<class ParamType>
PercentileAndMaximumFinder<ParamType>::~PercentileAndMaximumFinder()
{
    delete[] hists;
    delete[] cumSum;
    delete[] paramBounds;
}

template<class ParamType>
void PercentileAndMaximumFinder<ParamType>::processChainData(IttType start, IttType end, int sliceInterval)
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
            ++numPoints;
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
    if((!binned) ||  //the caller has not put in data
       (val<paramBounds[2*pIndex]) || //the value is below the range of parameters encountered
       (val>(paramBounds[2*pIndex]+paramBounds[2*pIndex+1]*(bCount+1)))) //the value is greater than the range of parameters encountered
    {
        return static_cast<ParamType>(-1);
    }
    //now translate the value back into a bin number
    int binNum = ((val - paramBounds[2*pIndex])/paramBounds[2*pIndex+1]);
    int cbinNum = (binNum + 1);
    //now linearly interpolate along the cumsum histogram to find the exact number of entries corresponding to that value
    //between the current location and the one before it to find the value
    //the value of the cell represents the value at the top edge of the cell,
    //the x values are then the high edge of the cell before and the current cell
    //the cell edge formula from paramBounds of binNum*binWidth+lowestBinEdge
    //gives the low edge for a bin, to get the high edge, add 1 to the bin number
    //however, the bin numbering is for hists, the cumSum arrays are offset up by
    //one bin, so subtract 1 from the bin number to get the high edge of the correct bin
    ParamType x1 = (((cbinNum-1)*paramBounds[2*pIndex+1])+paramBounds[2*pIndex]);
    ParamType x2 = ((cbinNum*paramBounds[2*pIndex+1])+paramBounds[2*pIndex]);
    ParamType y1 = cumSum[cbinNum-1];
    ParamType y2 = cumSum[cbinNum];
    //we have x1,y1 and x2,y2 so we can now find the slope and intercept of a straight-line connecting them
    ParamType m = ((y2-y1)/(x2-x1));
    ParamType b = (((x2*y1)-(x1*y2))/(x2-x1));
    //then the value we extract is the value the user wanted and simply comes from plugging in entries into the line formula
    ParamType entries = static_cast<int>((val*m)+b);
    //now we find the percentile by ratio
    return (static_cast<ParamType>(100)*(entries/static_cast<ParamType>(numPoints)));
}

template<class ParamType>
ParamType PercentileAndMaximumFinder<ParamType>::getValueFromPercentile(int pIndex, const ParamType& per)
{
    if((!binned) || (0>per) || (100<per)) return (paramBounds[2*pIndex] - 1e4);
    //first find out how many entries the given percentile is:
    int entries = (per/static_cast<ParamType>(100))*static_cast<ParamType>(numPoints);
    //get the cumulative sum array for the parameter
    int* cs = (cumSum + pIndex*cbCount);
    //do a binary search in the array to find the cell that contains the value (or the first cell whose value is above it)
    int first = 0;
    int last = (cbCount-1);
    if(entries == 0)
    {
        first = 1;
        last = 1;
    }
    else
    {
        while(first != last)
        {
            int mid = ((first + last)/2);
            if(cs[mid] >= entries)
            {//since cell 0 always contains 0 this condition cannot be true if mid = 0;
                if(cs[mid-1] <= entries)//if this is true we are at the correct cell
                {
                    first = mid;
                    last = mid;
                }
                else
                {//mid is too high a cell, try again
                    last = mid;
                }
            }
            else
            {// if we are here then mid was too low a cell, try again
                first = mid;
            }
        }
    }
    //now both first and last point to the correct location, do a linear interpolation
    //between the current location and the one before it to find the value
    //the value of the cell represents the value at the top edge of the cell,
    //the y values are then the high edge of the cell before and the current cell
    //the cell edge formula from paramBounds of binNum*binWidth+lowestBinEdge
    //gives the low edge for a bin, to get the high edge, add 1 to the bin number
    //however, the bin numbering is for hists, the cumSum arrays are offset up by
    //one bin, so subtract 1 from the bin number to get the high edge of the correct bin
    ParamType x1 = cs[last-1];
    ParamType x2 = cs[last];
    ParamType y1 = (((last-1)*paramBounds[2*pIndex+1])+paramBounds[2*pIndex]);
    ParamType y2 = ((last*paramBounds[2*pIndex+1])+paramBounds[2*pIndex]);
    //we have x1,y1 and x2,y2 so we can now find the slope and intercept of a straight-line connecting them
    ParamType m = ((y2-y1)/(x2-x1));
    ParamType b = (((x2*y1)-(x1*y2))/(x2-x1));
    //then the value we extract is the value the user wanted and simply comes from plugging in entries into the line formula
    return ((entries*m)+b);
}

template<class ParamType>
ParamType PercentileAndMaximumFinder<ParamType>::getValueOfPeak(int pIndex)
{
    if(!binned) return (paramBounds[2*pIndex] - 1e4);
    //jump to the histogram for the given parameter
    int* hist = (hists + pIndex*bCount);
    //find the bin with the largest content
    int binMax = 0;
    int maxVal = -1;
    for(int i=0; i<bCount; ++i)
    {
        if(hist[i] > maxVal)
        {
            maxVal = hist[i];
            binMax = i;
        }
    }
    //now translate that bin number into a bin center
    return (((binMax+0.5)*paramBounds[2*pIndex+1])+paramBounds[2*pIndex]);
}

template<class ParamType>
void PercentileAndMaximumFinder<ParamType>::writeHistogramsInCsvFormat(const std::string& fileNameBase)
{
    for(int i=0; i<pCount; ++i)
    {
        write1dHistCsv(i, fileNameBase);
        write1dCumSumCsv(i, fileNameBase);
    }
}

template<class ParamType>
void PercentileAndMaximumFinder<ParamType>::write1dHistCsv(int histNum, const std::string& fileNameBase)
{
    std::ostringstream namer;
    namer<<fileNameBase<<"_p"<<histNum<<".csv";
    std::ofstream out;
    out.open(namer.str().c_str());
    out<<"# Lines starting with a '#' in the first column are ignored\n";
    out<<"# X-axis: nbins, first bin low edge, last bin high edge\n";
    //write the axis parameters
    out<<bCount<<", "<<paramBounds[2*histNum]<<", "<<(paramBounds[2*histNum] + (bCount*paramBounds[2*histNum+1]))<<"\n";
    int offset = histNum*bCount;
    out<<"# bin number, value\n";
    //write the bin data
    for(int i=0; i<bCount; ++i)
    {
        out<<i<<", "<<hists[offset+i]<<"\n";
    }
    out<<std::flush;
    out.close();
}

template<class ParamType>
void PercentileAndMaximumFinder<ParamType>::write1dCumSumCsv(int histNum, const std::string& fileNameBase)
{
    std::ostringstream namer;
    namer<<fileNameBase<<"_cs_p"<<histNum<<".csv";
    std::ofstream out;
    out.open(namer.str().c_str());
    out<<"# Lines starting with a '#' in the first column are ignored\n";
    out<<"# X-axis: nbins, first bin low edge, last bin high edge\n";
    //write the axis parameters
    out<<bCount<<", "<<(paramBounds[2*histNum]-paramBounds[2*histNum+1])<<", "<<(paramBounds[2*histNum] + (bCount*paramBounds[2*histNum+1]))<<"\n";
    int offset = histNum*cbCount;
    out<<"# bin number, value\n";
    //write the bin data
    for(int i=0; i<cbCount; ++i)
    {
        out<<i<<", "<<cumSum[offset+i]<<"\n";
    }
    out<<std::flush;
    out.close();
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
    //reset the point counter
    numPoints = 0;
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
    binned = false;
}

template<class ParamType>
int PercentileAndMaximumFinder<ParamType>::sign(ParamType val)
{
    return (ParamType(0) < val) - (val < ParamType(0));
}

}
}

#endif //MCMCPP_ANALYSIS_PERCENTILES_H
