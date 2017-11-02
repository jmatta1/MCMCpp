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
#include<cmath>
#include<random>
// includes from other libraries
// includes from MCMC
#include"../Chain/ChainStepIterator.h"

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
 */

namespace Detail
{
const unsigned long long PiNumerator =   3141592653589793239ULL;///<Stores the numerator of Pi in rational form
const unsigned long long PiDenomenator = 1000000000000000000ULL;///<Stores the denomenator of Pi in rational form
}

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
    AutoCorrCalc(int numParams, int numWalkers) : paramCount(numParams), walkerCount(numWalkers)
    {acorrTimeList = new ParamType[paramCount]; randomWalkerIndices = new int[walkerCount]; for(int i=0; i<paramCount; ++i) acorrTimeList[i] = 0.0;}
    
    ~AutoCorrCalc()
    {delete[] acorrTimeList; delete[] randomWalkerIndices; if(acovFuncAvgArray!=nullptr) delete[] acovFuncAvgArray; if(acovFuncArray!=nullptr) delete[] acovFuncArray;
    if(interFuncArray!=nullptr) delete[] interFuncArray;}
    /*!
     * \brief allAutoCorrTime Calculates the auto correlation time for each parameter using the full set of walkers
     * \param start The iterator pointing at the start of the time series
     * \param end The iterator pointing past the last element in the time series
     * \param numSamples The number of samples, per walker, in the chain
     * 
     * Warning, this *can* be slow. Using fft methods it will take time proportional to
     * p*w*n*log2[n] where n is the number of samples to be used in the chain (the smallest power of 2
     * that is greater than or equal to the actual chain size), p is the number of parameters,
     * and w is the number of walkers.
     * 
     * Autocorrelation times that were calculated can be extracted using the retrieveAutoCorrelationTime function
     */
    void allAutoCorrTime(const IttType& start, const IttType& end, int numSamples);
    
    /*!
     * \brief sampleParameterAutoCorrTimes Calculates the autocorrellation for a given parameter using a subset of the walkers
     * \param numSamples The number of samples, per walker, in the chain
     * \param paramNumber The parameter index to calculate the autocorrelation time for
     * \param walkersToSelect The number of walkers to calculate the autocorrelation time for
     * \param keepPreviousWalkers If true, keeps the selection of walkers the acor time was previously calculated for
     * \return The autocorrelation time calculated for that parameter with that walker sampling
     * 
     * This too *can* be slow, using fft methods it will have time s*n*log2[n] where n has the same
     * meaning as in allAutoCorrTime and s is the number of walkers to sample. If fast is not set then
     * the time is proportional to s*n*log2[n]
     */
    ParamType sampleParameterAutoCorrTimes(const IttType& start, const IttType& end, int numSamples, int paramNumber, int walkersToSelect, bool keepPreviousWalkers=false)
    {checkScratchSizes(numSamples); return sampleParamAutoCorrTimesInternal(start, end, numSamples, paramNumber, walkersToSelect, keepPreviousWalkers);}
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
    ParamType sampleParamAutoCorrTimesInternal(const IttType& start, const IttType& end, int numSamples, int paramNumber, int walkersToSelect, bool keepPreviousWalkers);
    void genWalkerIndexList(int walkersToSelect);
    void checkScratchSizes(int numSamples);
    void genWalkerAutoCovFunc(const IttType& start, const IttType& end, int walkerNum, int numSamples, int paramNumber, int walkersToSelect);
    void makeCenteredWalkerChain(const IttType& start, const IttType& end, int walkerNum, int numSamples, int paramNumber);
    int bitReverse(int input);
    void fft();
    void copyBitInverseMagnitudes();
    void ifft();
    void averageAutocovarianceFunctions(int walkersToSelect);
    
    //General bookkeeping parameters
    int paramCount; ///<stores the number of parameters per sample
    int walkerCount; ///<stores the number of walkers in the chain
    
    //storage array for precalculated autocorrelation times
    ParamType* acorrTimeList; ///<stores the list of computed autocorrelation times calculated by allAutoCorrTime
    
    //general use random number generator stuff
    std::mt19937_64 engine;///<The base random number generator engine that is used for selecting walkers randomly
    std::normal_distribution<ParamType> normDist;///<The adapter that gives normally distributed real numbers with mean 0 and variance 1
    
    //Parameters for calculation of autocorellations
    ParamType* acovFuncAvgArray = nullptr; ///<Stores the sum of the autocovariance functions as they are calculated for every walker in the chain
    int acovSize = 0; ///<Stores the size of the autocovariance function arrays
    std::complex<ParamType>* acovFuncArray = nullptr; ///<Stores the autocovariance function calculated for a given walker in the chain
    std::complex<ParamType>* interFuncArray = nullptr; ///<stores the inverse fft generated in the first step of calculating the autocovariance function
    int scratchSize = 0; ///<Size of the acovFuncArray and interFuncArray arrays
    int fftSize = 0; ///<Stores the size of the fft
    int logFftSize = 0; ///<Stores the Log2 of the scratchSize
    int* randomWalkerIndices = nullptr; ///<Stores the array of randomly chosen walker indices
    int windowScaling = 4; ///<Minimum number of autocorrelation times to be processed to consider the result correct
    //has enough digits that it will be truncated to whatever precision is necessary (assuming the options are 32, 64, 80, or even 128 bit floats)
    const ParamType Pi = (static_cast<ParamType>(Detail::PiNumerator)/static_cast<ParamType>(Detail::PiDenomenator));
};

template<class ParamType>
void AutoCorrCalc<ParamType>::allAutoCorrTime(const IttType& start, const IttType& end, int numSamples)
{
    //first do the checks on scratch size
    checkScratchSizes(numSamples);
    //do the first one seperately to force the setting of the walker indice array
    acorrTimeList[0] = sampleParamAutoCorrTimesInternal(start, end, numSamples, 0, walkerCount, false);
    for(int i=1; i<= paramCount; ++i)
    {//simply apply the more limited autocorrelation time calculator multiple times, storing the result
        acorrTimeList[i] = sampleParamAutoCorrTimesInternal(start, end, numSamples, i, walkerCount, true);
    }
}

template<class ParamType>
ParamType AutoCorrCalc<ParamType>::sampleParamAutoCorrTimesInternal(const IttType& start, const IttType& end, int numSamples, int paramNumber, int walkersToSelect, bool keepPreviousWalkers)
{
    //now select the set of walkers whose autocorrelation functions are to be averaged
    if(!keepPreviousWalkers)
    {
        genWalkerIndexList(walkersToSelect);
    }
    //clear the autocovariance function average array
    for(int i=0; i<fftSize; ++i)
    {
        acovFuncAvgArray[i] = static_cast<ParamType>(0);
    }
    //now calculate the autocovariance function of the appropriate parameter in the selected chains
    for(int i=0; i<walkersToSelect; ++i)
    {
        //generate the autocorrelation function for a single walker
        genWalkerAutoCovFunc(start, end, randomWalkerIndices[i], numSamples, paramNumber, walkersToSelect);
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
void AutoCorrCalc<ParamType>::checkScratchSizes(int numSamples)
{
    //first check how many points we are using
    logFftSize = static_cast<int>(std::ceil(std::log2(numSamples)));
    int tempFftSize = (0x1UL << logFftSize);
    //now make sure that the storage for the autocovariance function is large enough
    if(acovSize < numSamples)
    {
        acovSize = numSamples;
        if(acovFuncAvgArray != nullptr) delete[] acovFuncAvgArray;
        acovFuncAvgArray = new ParamType[acovSize];
    }
    //now make sure that the storage for the inverted fft is large enough
    if(fftSize < tempFftSize)
    {
        fftSize = tempFftSize;
        if(acovFuncArray != nullptr) delete[] acovFuncArray;
        acovFuncArray = new std::complex<ParamType>[fftSize];
        if(interFuncArray != nullptr) delete[] interFuncArray;
        interFuncArray = new std::complex<ParamType>[fftSize];
    }
}

template<class ParamType>
void AutoCorrCalc<ParamType>::genWalkerAutoCovFunc(const IttType& start, const IttType& end, int walkerNum, int numSamples, int paramNumber, int walkersToSelect)
{
    // First calculate the average of the chain and center it
    makeCenteredWalkerChain(start, end, walkerNum, numSamples, paramNumber);
    // take the fft of the centered, zero-extended, walker chain 
    fft();
    // transfer/bitreverse and square the fft
    copyBitInverseMagnitudes();
    // do the inverse fast fourier transform of the norms of the forward fft
    ifft();
    //normalize the ifft and add it to the averaged data
    averageAutocovarianceFunctions(walkersToSelect);
}

template<class ParamType>
void AutoCorrCalc<ParamType>::averageAutocovarianceFunctions(int walkersToSelect)
{
    //This incorporates the normalization of the IFFT (1/N), the averaging factor (1/walkersToSelect), *and* the division by the autocovariance at lag 0
    //The since acovFuncArray[0] contains covFunc[t=0]*N (because it has not been normalized from the inverse FFT)
    //I can ignore dividing by 1/FFT size since that is already included
    ParamType normVal = (static_cast<ParamType>(1)/(static_cast<ParamType>(walkersToSelect)*(acovFuncArray[0].real())));
    for(int i=0; i<acovSize; ++i)
    {
        acovFuncAvgArray += (normVal*(acovFuncArray[i].real()));
    }
}

template<class ParamType>
void AutoCorrCalc<ParamType>::ifft()
{//simple implementation of Cooley-Tukey in reverse (negate frequency in exponent of freqStep and apply a factor of 1/N to everything after the fact
    //however this IFFT does not do the 1/N factor (that will be taken into account later,
    //during the averaging process, because otherwise it is an unnecessary loop through the array
    //which can be deferred to be combined with the necessary loop through the array to do the averaging)
    for(int s=0; s<logFftSize; ++s)    
    {
        int m1 = 0x1<<s;
        int m2 = m1 >> 1;
        std::complex<ParamType> freqStep(std::polar(static_cast<ParamType>(1), Pi/static_cast<ParamType>(m2)));
        for(unsigned int k=0; k<fftSize; k+=m1)
        {
            std::complex<ParamType> baseFreq(1, 0);
            for(unsigned int j=0; j<m2; ++j)
            {
                std::complex<ParamType> temp1 = (baseFreq * acovFuncArray[k+j+m2]);
                std::complex<ParamType> temp2 = acovFuncArray[k+j];
                acovFuncArray[k+j] = (temp2 + temp1);
                acovFuncArray[k+j+m2] = (temp2 - temp1);
                baseFreq *= freqStep;
            }
        }
    }
}

template<class ParamType>
void AutoCorrCalc<ParamType>::copyBitInverseMagnitudes()
{
    for(int i=0; i<fftSize; ++i)
    {
        acovFuncArray[bitReverse(i)] = std::norm(interFuncArray[i]);
    }
}

template<class ParamType>
void AutoCorrCalc<ParamType>::fft()
{//simple implementation of Cooley-Tukey
    for(int s=0; s<logFftSize; ++s)    
    {
        int m1 = 0x1<<s;
        int m2 = m1 >> 1;
        std::complex<ParamType> freqStep(std::polar(static_cast<ParamType>(1), -Pi/static_cast<ParamType>(m2)));
        for(unsigned int i=0; i<fftSize; i+=m1)
        {
            std::complex<ParamType> baseFreq(1, 0);
            for(unsigned int j=0; j<m2; ++j)
            {
                std::complex<ParamType> temp1 = (baseFreq * interFuncArray[k+j+m2]);
                std::complex<ParamType> temp2 = interFuncArray[k+j];
                interFuncArray[k+j] = (temp2 + temp1);
                interFuncArray[k+j+m2] = (temp2 - temp1);
                baseFreq *= freqStep;
            }
        }
    }
}

//Since chains might be long and have weird values, use the Kahan summation to reduce the floating point error
//Kahan summation can be optimized away by overly aggressive compilers, so it might be needed to turn the optimization down a bit for this function
//bottom of the file has code to wrap around this function if necessary
template<class ParamType>
void AutoCorrCalc<ParamType>::makeCenteredWalkerChain(const IttType& start, const IttType& end, int walkerNum, int numSamples, int paramNumber)
{
    ParamType sum = 0.0;
    ParamType compensation = 0.0;
    int offset = (walkerNum*paramCount + paramNumber);
    int index = 0;
    for(IttType itt(start); itt != end; ++itt)
    {
        //Sample copying
        ParamType temp = (*itt)[offset];
        interFuncArray[bitReverse(index)] = temp;
        ++index;
        //Kahan addition for increased accuracy in calculating the average
        ParamType temp1 = (temp - compensation);
        ParamType temp2 = sum + temp1;  // If sum is big and temp1 is small low order digits of (*itt)[offset] are lost
        compensation = (temp2 - sum) - temp1; //restore the low order digits of (*itt)[offset] to try again with later
        sum = temp2;
    }
    ParamType avg = (sum/static_cast<ParamType>(numSamples));
    for(int i=0; i<numSamples; ++i)
    {
        interFuncArray -= avg;
    }
    //pretend that we have zero extended the data array before we did this whole copy, subtract, and bit-reversal thing
    for(; index<fftSize; ++index)
    {
        interFuncArray[bitReverse(index)] = 0.0;
    }
}

template<class ParamType>
int AutoCorrCalc<ParamType>::bitReverse(int input)
{
    int output = 0;
    for(int i=0; i<logFftSize; ++i)
    {
        output <<= 1;
        output |= (input & 0x1);
        input >>= 1;
    }
    return output;
}

template<class ParamType>
void AutoCorrCalc<ParamType>::genWalkerIndexList(int walkersToSelect)
{
    if(walkerCount == walkersToSelect)
    {
        for(int i=0; i<=walkerCount; ++i)
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
/*
 * Code to wrap around function to reduce optimization level
 * useful for functions where aggressive transitivity rules might collapse the whole scheme,
 * like for Kahan summation
 * 
#ifdef __GNUC__
#pragma GCC push_options
#pragma GCC optimize ("O2")
#endif
#ifdef __INTEL_COMPILER
#pragma optimize off
#endif
#ifdef __clang__
#pragma clang optimize off
#endif
  //Put function here
#ifdef __GNUC__
#pragma GCC pop_options
#endif 
#ifdef __INTEL_COMPILER
#pragma optimize on
#endif
#ifdef __clang__
#pragma clang optimize on
#endif
*/
}
}
#endif  //MCMC_ANALYSIS_AUTOCORRCALC_H
