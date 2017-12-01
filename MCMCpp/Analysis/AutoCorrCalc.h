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
#ifndef MCMC_ANALYSIS_AUTOCORRCALC_H
#define MCMC_ANALYSIS_AUTOCORRCALC_H
// includes for C system headers
// includes for C++ system headers
#include<complex>//needed for the FFT and iFFT
#include<cmath>//needed for ceiling and log 2
#include<algorithm>//used for
#include<random>//needed for normal distribution
// includes from other libraries
#include"../Utility/pcg-cpp/include/pcg_random.hpp"
// includes from MCMC
#include"../Chain/ChainStepIterator.h"
#include"Detail/AutoCov.h"

namespace MCMC
{
namespace Analysis
{

/*!
 * @class AutoCorrCalc
 * @ingroup Analysis
 * @brief A class to calculate the autocorrelation times of the parameters that were MCMCed
 * @author James Till Matta
 * 
 * @tparam ParamType The floating point type for which this calculation will be carried out
 */
template<class ParamType>
class AutoCorrCalc
{
public:
    typedef Chain::ChainStepIterator<ParamType> IttType;
    /*!
     * \brief AutoCorrCalc constructs a new AutoCorrCalc object
     * \param numParams The number of parameters in each sample
     * \param numWalkers The number of walkers in the ensemble
     */
    AutoCorrCalc(int numParams, int numWalkers) :
        engine(std::random_device()()), normDist(static_cast<ParamType>(0),static_cast<ParamType>(1)),
        autoCovCalc(), paramCount(numParams), walkerCount(numWalkers)
    {   acorrTimeList = new ParamType[paramCount]; randomWalkerIndices = new int[walkerCount];
        chainAverages = new ParamType[paramCount*numWalkers];  std::fill_n(acorrTimeList, paramCount, static_cast<ParamType>(0));}
    
    ~AutoCorrCalc()
    {delete[] acorrTimeList; delete[] randomWalkerIndices; delete[] chainAverages; if(acovFuncAvgArray!=nullptr) delete[] acovFuncAvgArray;
        if(acovFuncArray!=nullptr) delete[] acovFuncArray;}

    /*!
     * \brief Deleted copy constructor
     */
    AutoCorrCalc(const AutoCorrCalc<ParamType>& rhs) = delete;
    
    /*!
     * \brief Deleted assignment operator
     */
    AutoCorrCalc<ParamType>& operator=(const AutoCorrCalc<ParamType>& rhs) = delete;
    
    /*!
     * \brief calcAutoCorrTimes Calculates the auto correlation time for each parameter using the full set of walkers
     * \param start The iterator pointing at the start of the time series
     * \param end The iterator pointing past the last element in the time series
     * \param numSamples The number of samples, per walker, in the chain
     * \param numWalkersToUse The number of randomly selected walkers (from the full ensemble) to use
     * 
     * Warning, this *can* be slow. Using fft methods it will take time proportional to
     * p*w*n*log2[n] where n is the number of samples to be used in the chain (the smallest power of 2
     * that is greater than or equal to the actual chain size), p is the number of parameters,
     * and w is the number of walkers.
     * 
     * Autocorrelation times that were calculated can be extracted using the retrieveAutoCorrelationTime function
     */
    void calcAutoCorrTimes(const IttType& start, const IttType& end, int numSamples, int numWalkersToUse = 0);

    /*!
     * \brief setAutoCorrParameters Sets the parameters for calculation of the autocorrelation time
     * \param scaleFactor The minimum required number of autocorrelation times the algorithm needs to examine
     */
    void setAutoCorrScaleFactor(int scaleFactor=5){windowScaling = scaleFactor;}
    /*!
     * \brief getAutoCorrelationTime retrieves the calculated autocorrelation time for parameter number paramIndex
     * \param paramIndex The index of the parameter [0, numParameter)
     * \return The autocorrelation time in samples for parameter # paramIndex
     */
    ParamType retrieveAutoCorrelationTime(int paramIndex){return acorrTimeList[paramIndex];}
    
private:
    //"Behind the scenes" functions that do some menial work and heavy lifting
    ParamType sampleParamAutoCorrTimes(const IttType& start, const IttType& end, int numSamples, int paramNumber, int walkersToSelect);
    void genWalkerIndexList(int walkersToSelect);
    void checkScratchSizes(int numSamples);
    void calculateChainAverages(const IttType& start, const IttType& end, int numSamples, int numWalkersToUse);
    void transferWalker(const IttType& start, const IttType& end, int paramNumber, int walkerNumber);
    void averageAutocovarianceFunctions(int walkersToSelect);
    
    //general use random number generator stuff
    pcg32 engine;///<The base random number generator engine that is used for selecting walkers randomly
    std::normal_distribution<ParamType> normDist;///<The adapter that gives normally distributed real numbers with mean 0 and variance 1
    
    Detail::AutoCov<ParamType> autoCovCalc; ///<Performs the calculation of the autocovariance function
    
    //storage array for precalculated autocorrelation times
    ParamType* acorrTimeList = nullptr; ///<stores the list of computed autocorrelation times calculated by allAutoCorrTime
    
    //Parameters for calculation of autocorellations
    ParamType* acovFuncAvgArray = nullptr; ///<Stores the sum of the autocovariance functions as they are calculated for every walker in the chain
    ParamType* acovFuncArray = nullptr; ///<Stores the chain, and the autocovariance function that is calculated by the object
    ParamType* chainAverages = nullptr; ///<Stores the computed averages of the chains selected for calculation
    int* randomWalkerIndices = nullptr; ///<Stores the array of randomly chosen walker indices
    int acovSize = 0; ///<Stores the size of the autocovariance function arrays
    int windowScaling = 4; ///<Minimum number of autocorrelation times to be processed to consider the result correct

    //General bookkeeping parameters
    int paramCount; ///<stores the number of parameters per sample
    int walkerCount; ///<stores the number of walkers in the chain
};

template<class ParamType>
void AutoCorrCalc<ParamType>::calcAutoCorrTimes(const IttType& start, const IttType& end, int numSamples, int numWalkersToUse)
{
    if(numWalkersToUse == 0) numWalkersToUse = walkerCount;
    //first do the checks on scratch size
    checkScratchSizes(numSamples);
    //now generate the list of walkers
    genWalkerIndexList(numWalkersToUse);
    //now calculate the chain averages
    calculateChainAverages(start, end, numSamples, numWalkersToUse);
    //do the first one seperately to force the setting of the walker indice array
    acorrTimeList[0] = sampleParamAutoCorrTimes(start, end, numSamples, 0, numWalkersToUse);
    //simply apply the more limited single parameter autocorrelation time calculator multiple times, storing the result
    for(int i=1; i< paramCount; ++i)
    {
        acorrTimeList[i] = sampleParamAutoCorrTimes(start, end, numSamples, i, numWalkersToUse);
    }
}

template<class ParamType>
ParamType AutoCorrCalc<ParamType>::sampleParamAutoCorrTimes(const IttType& start, const IttType& end, int numSamples, int paramNumber, int walkersToSelect)
{
    //clear the autocovariance function average array
    std::fill_n(acovFuncAvgArray, acovSize, static_cast<ParamType>(0));
    
    //now calculate the autocovariance function for every selected walker and add it to the average
    for(int i=0; i<walkersToSelect; ++i)
    {
        transferWalker(start, end, paramNumber, randomWalkerIndices[i]);
        autoCovCalc.calcNormAutoCov(acovFuncArray,chainAverages[randomWalkerIndices[i]*paramCount + paramNumber], numSamples);
        averageAutocovarianceFunctions(walkersToSelect);
    }
    //now that we have the averaged autocovariance function, extract the autocorrelation time from the cumulative sum
    ParamType autoCorrSum = acovFuncAvgArray[0]; // initialize to this so we only count the first cell once
    ParamType factor = static_cast<ParamType>(windowScaling);
    //because of definitions, the autocorrelation time cannot possibly be less than 1 since acovFuncAvgArray[0] == 1
    for(int i=1; i<acovSize; ++i)
    {
        //Add the next term to the autocovariance function
        autoCorrSum += (static_cast<ParamType>(2)*acovFuncAvgArray[i]);
        //check if our window size surpasses the currently estimated autocorrelation time
        //if it does return the currently estimated autocorrelation time
        if(i > factor*autoCorrSum)
        {
            return autoCorrSum;
        }
        //otherwise keep looping
    }
    // if we got to here then the autocorrelation time never converged for a
    // given safety factor, return the *negative* of the final value in the summation
    return -autoCorrSum;
}

template<class ParamType>
void AutoCorrCalc<ParamType>::averageAutocovarianceFunctions(int walkersToSelect)
{
    ParamType normVal = (static_cast<ParamType>(1)/(static_cast<ParamType>(walkersToSelect)));
    std::transform(acovFuncArray, acovFuncArray+acovSize, acovFuncAvgArray, acovFuncAvgArray,
                   [&normVal] (const ParamType& acv, const ParamType& avg) -> ParamType {return (avg+(normVal*acv));});
}

template<class ParamType>
void AutoCorrCalc<ParamType>::transferWalker(const IttType& start, const IttType& end, int paramNumber, int walkerNumber)
{
    int offset = (walkerNumber*paramCount+paramNumber);
    IttType itt(start);
    std::transform(itt, end, acovFuncArray, acovFuncArray,
                   [&offset](ParamType* ptr, const ParamType& acv) -> ParamType {return (ptr[offset]+acv);});
}

template<class ParamType>
void AutoCorrCalc<ParamType>::calculateChainAverages(const IttType& start, const IttType& end, int numSamples, int numWalkersToUse)
{
    ParamType norm = (static_cast<ParamType>(1)/static_cast<ParamType>(numSamples));
    for(IttType itt(start); itt != end; ++itt)
    {
        for(int i=0; i<numWalkersToUse; ++i)
        {
            int offset = randomWalkerIndices[i]*paramCount;
            int endset = (offset+paramCount);
            for(int j=offset; j<endset; ++j)
            {
                chainAverages[j] += (norm*((*itt)[j]));
            }
        }
    }
}

template<class ParamType>
void AutoCorrCalc<ParamType>::genWalkerIndexList(int walkersToSelect)
{
    if(walkerCount == walkersToSelect)
    {
        for(int i=0; i<walkerCount; ++i)
        {
            randomWalkerIndices[i] = i;
        }
    }
    else
    {
        int numSelected = 0;
        int totalInputExamined = 0;
        ParamType randUniform;
        while(numSelected < walkersToSelect)
        {
            randUniform = normDist(engine);
            if( ((walkerCount-totalInputExamined)*randUniform) >= (walkerCount-numSelected))
            {
                ++totalInputExamined;
            }
            else
            {
                randomWalkerIndices[numSelected] = totalInputExamined;
                ++totalInputExamined;
                ++numSelected;
            }
        }
    }
}

template<class ParamType>
void AutoCorrCalc<ParamType>::checkScratchSizes(int numSamples)
{
    //now make sure that the storage for the autocovariance function is large enough
    if(acovSize < numSamples)
    {
        acovSize = numSamples;
        if(acovFuncAvgArray != nullptr) delete[] acovFuncAvgArray;
        acovFuncAvgArray = new ParamType[acovSize];
        if(acovFuncArray != nullptr) delete[] acovFuncArray;
        acovFuncArray = new ParamType[acovSize];
    }
}

}
}
#endif  //MCMC_ANALYSIS_AUTOCORRCALC_H
