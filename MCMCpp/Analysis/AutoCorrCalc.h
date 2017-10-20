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
#include<complex>//needed for the DFT and FFT
#include<cmath>
// includes from other libraries
// includes from MCMC

namespace MarkovChainMonteCarlo
{
namespace Analysis
{
/*!
 * @class AutoCorrCalc
 * @ingroup Utility
 * @brief A class to calculate the autocorrelation times of the parameters that were MCMCed
 * @author James Till Matta
 * 
 * @tparam ParamType The floating point type for which this calculation will be carried out
 * @tparam IttType The type of iterator for the time series that it will be working with
 */
template<class ParamType, class IttType>
class AutoCorrCalc
{
    
    /*!
     * \brief AutoCorrCalc constructs a new AutoCorrCalc object
     * \param numParams The number of parameters in each sample
     * \param numWalkers The number of walkers in the ensemble
     */
    AutoCorrCalc(int numParams, int numWalkers) : paramCount(numParams), walkerCount(numWalkers)
    {acorrTimeList = new ParamType[paramCount]; randomWalkerIndices = new int[walkerCount]; for(int i=0; i<paramCount; ++i) acorrTimeList[i] = 0.0;}
    
    ~AutoCorrCalc()
    {delete[] acorrTimeList; delete[] randomWalkerIndices; if(acovFuncAvgArray!=nullptr) delete[] acovFuncAvgArray; if(acovFuncArray!=nullptr) delete[] acovFuncArray;
    if(interFuncArray!=nullptr) delete[] interFuncArray;}
    /*!
     * \brief allAutoCorrTime Calculates the auto correlation time for each parameter using the full set of walkers
     * \param numSamples The number of samples, per walker, in the chain
     * 
     * Warning, this *can* be slow. Using fft methods it will take time proportional to
     * p*w*n*log2[n] where n is the number of samples to be used in the chain (the smallest power of 2
     * that is greater than or equal to the actual chain size), p is the number of parameters,
     * and w is the number of walkers.
     * 
     * Autocorrelation times that were calculated can be extracted using the retrieveAutoCorrelationTime function
     */
    void allAutoCorrTime(const MarkovChainMonteCarlo::Utility::IttType& start, const MarkovChainMonteCarlo::Utility::IttType& end, int numSamples);
    
    /*!
     * \brief sampleParameterAutoCorrTimes Calculates the autocorrellation for a given parameter using a subset of the walkers
     * \param numSamples The number of samples, per walker, in the chain
     * \param paramNumber The parameter index to calculate the autocorrelation time for
     * \param numWalkers The number of walkers to calculate the autocorrelation time for
     * \param randomizeWalkers If true, randomly selects walkersToSample walkers instead of evenly dividing the walkers
     * \return The autocorrelation time calculated for that parameter with that walker sampling
     * 
     * This too *can* be slow, using fft methods it will have time s*n*log2[n] where n has the same
     * meaning as in allAutoCorrTime and s is the number of walkers to sample. If fast is not set then
     * the time is proportional to s*n*log2[n]
     */
    ParamType sampleParameterAutoCorrTimes(const IttType& start, const IttType& end, int numSamples, int paramNumber, int numWalkers, bool randomizeWalkers=false);
    /*!
     * \brief setAutoCorrParameters Sets the parameters for calculation of the autocorrelation time
     * \param minAutoCorrTimes The minimum required number of autocorrelation times the algorithm needs to examine
     * \param step The increase in window size for each iteration of the algorithm
     * \param loWin The minimum window size
     * \param hiWin The maximum window size
     * \param fast If true, only use the first power of two samples to accellerate calculation using FFT
     */
    void setAutoCorrParameters(int minAutoCorrTimes=10, int step=1, int loWin=10, int hiWin=10000);
    /*!
     * \brief getAutoCorrelationTime retrieves the calculated autocorrelation time for parameter number paramIndex
     * \param paramIndex The index of the parameter [0, numParameter)
     * \return The autocorrelation time in samples for parameter # paramIndex
     */
    ParamType retrieveAutoCorrelationTime(int paramIndex){return acorrTimeList[paramIndex];}
    
private:
    ParamType* acorrTimeList; ///<stores the list of computed autocorrelation times calculated by allAutoCorrTime
    int paramCount; ///<stores the number of parameters per sample
    int walkerCount; ///<stores the number of walkers in the chain
    ParamType* acovFuncAvgArray = nullptr; ///<Stores the sum of the autocovariance functions as they are calculated for every walker in the chain
    ParamType* acovFuncArray = nullptr; ///<Stores the autocovariance function calculated for a given walker in the chain
    int acovSize = 0; ///<Stores the size of the autocovariance function arrays
    std::complex<ParamType>* interFuncArray = nullptr; ///<stores the inverse fft generated in the first step of calculating the autocovariance function
    int scratchSize = 0; ///<Size of the acovFuncSumArray, acovFuncArray, and interFuncArray arrays
    int* randomWalkerIndices = nullptr; ///<Stores the array of randomly chosen walker indices
    int minAcorTimes = 10; ///<Minimum number of autocorrelation times to be processed to consider the result correct
    int winStepSize = 1; ///<Size of step for increasing window size
    int minWinSize = 10; ///<minimum window size for the algorithm
    int maxWinSize = 10000; ///<maximum window size for the algorithm
};

template<class ParamType, class IttType>
void AutoCorrCalc<ParamType, IttType>::allAutoCorrTime(const IttType& start, const IttType& end, int numSamples)
{
    for(int i=0; i<= paramCount; ++i)
    {//simply apply the more limited autocorrelation time calculator multiple times, storing the result
        acorrTimeList[i] = sampleParameterAutoCorrTimes(start, end, numSamples, i, walkerCount);
    }
}

template<class ParamType, class IttType>
void AutoCorrCalc<ParamType, IttType>::setAutoCorrParameters(int minAutoCorrTimes=10, int step=1, int loWin=10, int hiWin=10000)
{
    minAcorTimes = minAutoCorrTimes;
    winStepSize = step;
    minWinSize = loWin;
    maxWinSize = hiWin;
}


template<class ParamType, class IttType>
ParamType AutoCorrCalc<ParamType, IttType>::
sampleParameterAutoCorrTimes(const IttType& start, const IttType& end, int numSamples, int paramNumber, int numWalkers, bool randomizeWalkers=false)
{
    //first check how many points we are using
    int fftSize = (0x01 << static_cast<int>(std::ceiling(std::log2(numSamples))));
    //now make sure that the storage for the autocovariance function is large enough
    if(acovSize < numSamples)
    {
        acovSize = numSamples;
        if(acovFuncArray != nullptr) delete[] acovFuncArray;
        if(acovFuncAvgArray != nullptr) delete[] acovFuncAvgArray;
        acovFuncArray = new ParamType[acovSize];
        acovFuncAvgArray = new ParamType[acovSize];
    }
    //now make sure that the storage for the inverted fft is large enough
    if(scratchSize < fftSize)
    {
        scratchSize = fftSize;
        if(interFuncArray != nullptr) delete[] interFuncArray;
        interFuncArray = new ParamType[scratchSize];
    }
    //now select the set of walkers whose autocorrelation functions are to be averaged
    int numSelected = 0;
    int totalInputExamined = 0;
    ParamType randUniform;
    while(numSelected < numPoints)
    {
        randUniform = prng.getUniformReal();
        if( ((numWalkers-totalInputExamined)*randUniform) >= (numWalkers-numSelected))
        {
            ++totalInputExamined;
        }
        else
        {
            walkerIndices[numSelected] = totalInputExamined;
            ++totalInputExamined;
            ++numSelected;
        }
    }
}

}
}
#endif  //MCMC_ANALYSIS_AUTOCORRCALC_H
