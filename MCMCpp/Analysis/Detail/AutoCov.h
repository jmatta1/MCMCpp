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
#ifndef MCMC_ANALYSIS_DETAIL_AUTOCOV_H
#define MCMC_ANALYSIS_DETAIL_AUTOCOV_H
// includes for C system headers
// includes for C++ system headers
#include<stdlib.h>//needef for aligned allocation
#include<complex>//needed for the FFT and iFFT
#include<cmath>//needed for ceiling and log 2
// includes from other libraries
// includes from MCMC

namespace MarkovChainMonteCarlo
{
namespace Analysis
{

namespace Detail
{
const unsigned long long PiNumerator =   3141592653589793239ULL;///<Stores the numerator of Pi in rational form
const unsigned long long PiDenomenator = 1000000000000000000ULL;///<Stores the denomenator of Pi in rational form
const unsigned int AlignmentLength = 64;///<Stores the memory boundary to force memory alignment to, 64 is sufficient for cache lines and up to the 256-bit AVX instructions, 128 will handle AVX-512 instructions as well
/*!
 * @class AutoCov
 * @ingroup AnalysisDetail
 * @brief A class to calculate the autocovariance function using fft from the centered time series of a parameter for a walker
 * @author James Till Matta
 * 
 * @tparam ParamType The floating point type for which this calculation will be carried out
 */
template<class ParamType>
class AutoCov
{
public:
    /*!
     * \brief AutoCov Construct an uninitialized autocovariance calculator object, initialization occurs at the first calculation carried out, reinitialization occurs
     */
    AutoCov(){}
    
    ~AutoCov();
    
    /*!
     * \brief calcNormAutoCov Calculates the normalized autocovariance for a centered chain and places it in the input
     * \param centeredChain Pointer to the start of an array containing the chain which will also contain the output when finished
     * \param avg Average value of the chain
     * \param chainLength Length of the chain in samples
     */
    void calcNormAutoCov(ParamType* chain, const ParamType& avg, int chainLength);
    
private:
    
    /*!
     * \brief normalizeAndCopyAutoCov Transfers the normalized autocovariance function into the given array
     * \param output the array to hold the function
     * \param chainLength the length of the original chain
     */
    void normalizeAndCopyAutoCov(ParamType* output, int chainLength);
    
    /*!
     * \brief performInverseFft Calculates the inverse fft of the power spectrum of the chain, neglects the 1/N normalization factor
     */
    void performInverseFft();
    
    /*!
     * \brief makeBitReversedPowerSpectrum Takes the result of the for forward fft and copies its magnitude^2 in a bit reversed index fashion into the secondTransform array
     */
    void copyBitReversedPowerSpectrum();
    
    /*!
     * \brief performForwardFft Calculates the forward FFT of the chain to extract the chains amplitude spectrum
     */
    void performForwardFft();
    
    /*!
     * \brief copyCenterAndBitReverseChain Copies the chain, in bitreversed order to firstTransform while subtracting the average and zero extending the chain
     * \param chain The chain to be copied
     * \param avg The average to be subtracted
     * \param chainLength The length of the chain in question
     */
    void copyCenterAndBitReverseChain(const ParamType* chain, const ParamType& avg, int chainLength);
    
    /*!
     * \brief bitReverse Reverses the bits of an index for use with the iterative bit reversed form of Cooley-Tukey
     * \param input The index in question
     * \return The bit reversed integer
     */
    inline unsigned int bitReverse(unsigned int input) const;
    
    //Initialization / Reinitialization functions
    /*!
     * \brief checkSizeAndHandleChanges Takes the size of the chain and determines if arrays are the correct size and reallocate them if they are not
     * \param chainLength The length of the chain that needs to be accomodated
     */
    void checkSizeAndHandleChanges(int chainLength);
    
    /*!
     * \brief calculateRootsOfUnity generates the roots of unity that are needed as "twiddle factors"
     */
    void calculateRootsOfUnity();
    
    /*!
     * \brief findNextPowerOfTwo Takes a size and finds the first power of two that would be greater than it
     * \param chainLength The length of the chain to "upsize"
     * \return The log2 of the next larger power of two
     */
    inline int findNextPowerOfTwo(int chainLength) const;
    
    int fftSize = 0; ///<stores the number of points in the FFT
    int lgFftSize = -1; ///< stores the log2 of the number of points in the fft size
    ParamType* rootsOfUnityReal = nullptr; ///<Stores the real part of e^((+/-)i*2*n/N) twiddle factors used in fft
    ParamType* negRootsOfUnityImag = nullptr; ///<Stores the imaginary part of e^(i*2*pi*k*n/N) twiddle factors used in ifft
    ParamType* posRootsOfUnityImag = nullptr; ///<Stores the imaginary part of e^(i*2*pi*k*n/N) twiddle factors used in ifft
    
    ParamType* firstTransformReal = nullptr; ///<Where the real part of the chain is initially copied to undergo its first transform (and for its conversion to a power spectrum)
    ParamType* firstTransformImag = nullptr; ///<Where the imaginary part of the chain is initially copied to undergo its first transform (and for its conversion to a power spectrum)
    ParamType* secondTransformReal = nullptr; ///<Where the real part of the chain undergoes its second transform
    ParamType* secondTransformImag = nullptr; ///<Where the imaginary part of the chain undergoes its second transform
    
    const ParamType TwoPi = (static_cast<ParamType>(2)*(static_cast<ParamType>(PiNumerator)/static_cast<ParamType>(PiDenomenator)));///<Storage of 2*pi
};

template<class ParamType>
AutoCov<ParamType>::~AutoCov()
{
    if(rootsOfUnityReal != nullptr) free(rootsOfUnityReal);
    if(negRootsOfUnityImag != nullptr) free(negRootsOfUnityImag);
    if(posRootsOfUnityImag != nullptr) free(posRootsOfUnityImag);
    if(firstTransformReal != nullptr) free(firstTransformReal);
    if(firstTransformImag != nullptr) free(firstTransformImag);
    if(secondTransformReal != nullptr) free(secondTransformReal);
    if(secondTransformImag != nullptr) free(secondTransformImag);
}

template<class ParamType>
void AutoCov<ParamType>::calcNormAutoCov(ParamType* chain, const ParamType& avg, int chainLength)
{
    checkSizeAndHandleChanges(chainLength);
    copyCenterAndBitReverseChain(chain, avg, chainLength);
    performForwardFft();
    copyBitReversedPowerSpectrum();
    performInverseFft();
    normalizeAndCopyAutoCov(chain, chainLength);
}

template<class ParamType>
void AutoCov<ParamType>::normalizeAndCopyAutoCov(ParamType* output, int chainLength)
{
    ParamType normValue = (static_cast<ParamType>(1)/secondTransformReal[0]);
    for(int i=0; i<chainLength; ++i)
    {
        output[i] = (normValue * secondTransformReal[i]);
    }
}

template<class ParamType>
void AutoCov<ParamType>::performInverseFft()
{
    alignas(16) ParamType tDat1[2] = {0.0, 0.0};
    alignas(16) ParamType tDat2[2] = {0.0, 0.0};
    for(int s=1; s<=lgFftSize; ++s)
    {
        int m1 = 0x1<<s;
        int m2 = m1 >> 1;
        int twiddleStride = (fftSize>>s);
        for(unsigned int k=0; k<fftSize; k+=m1)
        {
            int twiddleIndex = 0;
            for(unsigned int j=0; j<m2; ++j)
            {
                tDat1[0] = rootsOfUnityReal[twiddleIndex]*firstTransformReal[k+j+m2] - posRootsOfUnityImag[twiddleIndex]*firstTransformImag[k+j+m2];
                tDat1[1] = rootsOfUnityReal[twiddleIndex]*firstTransformImag[k+j+m2] + posRootsOfUnityImag[twiddleIndex]*firstTransformReal[k+j+m2];
                //std::complex<ParamType> temp1 = (negRootsOfUnity[twiddleIndex] * firstTransform[k+j+m2]);
                tDat2[0] = firstTransformReal[k+j];
                tDat2[1] = firstTransformImag[k+j];
                //std::complex<ParamType> temp2 = firstTransform[k+j];
                firstTransformReal[k+j] = (tDat2[0] + tDat1[0]);
                firstTransformImag[k+j] = (tDat2[1] + tDat1[1]);
                //firstTransform[k+j] = (temp2 + temp1);
                firstTransformReal[k+j+m2] = (tDat2[0] - tDat1[0]);
                firstTransformImag[k+j+m2] = (tDat2[1] - tDat1[1]);
                //firstTransform[k+j+m2] = (temp2 - temp1);
                twiddleIndex += twiddleStride;
            }
        }
    }
}

template<class ParamType>
void AutoCov<ParamType>::copyBitReversedPowerSpectrum()
{
    for(int i=0; i<fftSize; ++i)
    {
        int bitRev = bitReverse(i);
        secondTransformReal[bitRev] = ((firstTransformReal[i]*firstTransformReal[i]) +
                                       (firstTransformImag[i]*firstTransformImag[i]));
        secondTransformReal[bitRev] = static_cast<ParamType>(0);
    }
}

template<class ParamType>
void AutoCov<ParamType>::performForwardFft()
{//implementation of Cooley-Tukey
    alignas(16) ParamType tDat1[2] = {0.0, 0.0};
    alignas(16) ParamType tDat2[2] = {0.0, 0.0};
    for(int s=1; s<=lgFftSize; ++s)
    {
        int m1 = 0x1<<s;
        int m2 = m1 >> 1;
        int twiddleStride = (fftSize>>s);
        for(unsigned int k=0; k<fftSize; k+=m1)
        {
            int twiddleIndex = 0;
            for(unsigned int j=0; j<m2; ++j)
            {
                tDat1[0] = rootsOfUnityReal[twiddleIndex]*firstTransformReal[k+j+m2] - negRootsOfUnityImag[twiddleIndex]*firstTransformImag[k+j+m2];
                tDat1[1] = rootsOfUnityReal[twiddleIndex]*firstTransformImag[k+j+m2] + negRootsOfUnityImag[twiddleIndex]*firstTransformReal[k+j+m2];
                //std::complex<ParamType> temp1 = (negRootsOfUnity[twiddleIndex] * firstTransform[k+j+m2]);
                tDat2[0] = firstTransformReal[k+j];
                tDat2[1] = firstTransformImag[k+j];
                //std::complex<ParamType> temp2 = firstTransform[k+j];
                firstTransformReal[k+j] = (tDat2[0] + tDat1[0]);
                firstTransformImag[k+j] = (tDat2[1] + tDat1[1]);
                //firstTransform[k+j] = (temp2 + temp1);
                firstTransformReal[k+j+m2] = (tDat2[0] - tDat1[0]);
                firstTransformImag[k+j+m2] = (tDat2[1] - tDat1[1]);
                //firstTransform[k+j+m2] = (temp2 - temp1);
                twiddleIndex += twiddleStride;
            }
        }
    }
}

template<class ParamType>
void AutoCov<ParamType>::copyCenterAndBitReverseChain(const ParamType* chain, const ParamType& avg, int chainLength)
{
    for(int i=0; i<chainLength; ++i)
    {
        int bitRev = bitReverse(i);
        firstTransformReal[bitRev] = (chain[i] - avg);
        firstTransformImag[bitRev] = static_cast<ParamType>(0);
    }
    
    for(int i=chainLength; i<fftSize; ++i)
    {
        int bitRev = bitReverse(i);
        firstTransformReal[bitRev] = static_cast<ParamType>(0);
        firstTransformImag[bitRev] = static_cast<ParamType>(0);
    }
}

template<class ParamType>
unsigned int AutoCov<ParamType>::bitReverse(unsigned int input) const
{
    //swap even and odd bits
    input = ((input >> 1) & 0x55555555) | ((input & 0x55555555) << 1);
    input = ((input >> 2) & 0x33333333) | ((input & 0x33333333) << 2);
    input = ((input >> 4) & 0x0F0F0F0F) | ((input & 0x0F0F0F0F) << 4);
    input = ((input >> 8) & 0x00FF00FF) | ((input & 0x00FF00FF) << 8);
    input = (input >> 16) | (input << 16);
    //now we have the 32-bit number reversed, bit shift it down by the number of unused bits
    input >>= (32 - lgFftSize);
    return input;
}

template<class ParamType>
void AutoCov<ParamType>::checkSizeAndHandleChanges(int chainLength)
{
    int tempLgFftSize = findNextPowerOfTwo(chainLength);
    if(tempLgFftSize != lgFftSize)
    {
        lgFftSize = tempLgFftSize;
        fftSize = (1 << lgFftSize);
        if(rootsOfUnityReal    != nullptr) free(rootsOfUnityReal);
        if(negRootsOfUnityImag != nullptr) free(negRootsOfUnityImag);
        if(posRootsOfUnityImag != nullptr) free(posRootsOfUnityImag);
        if(firstTransformReal  != nullptr) free(firstTransformReal);
        if(firstTransformImag  != nullptr) free(firstTransformImag);
        if(secondTransformReal != nullptr) free(secondTransformReal);
        if(secondTransformImag != nullptr) free(secondTransformImag);
        int twiddleSize = (fftSize/2);
        rootsOfUnityReal    = reinterpret_cast<ParamType*>(aligned_alloc(AlignmentLength, sizeof(ParamType)*twiddleSize));
        negRootsOfUnityImag = reinterpret_cast<ParamType*>(aligned_alloc(AlignmentLength, sizeof(ParamType)*twiddleSize));
        posRootsOfUnityImag = reinterpret_cast<ParamType*>(aligned_alloc(AlignmentLength, sizeof(ParamType)*twiddleSize));
        firstTransformReal  = reinterpret_cast<ParamType*>(aligned_alloc(AlignmentLength, sizeof(ParamType)*fftSize));
        firstTransformImag  = reinterpret_cast<ParamType*>(aligned_alloc(AlignmentLength, sizeof(ParamType)*fftSize));
        secondTransformReal = reinterpret_cast<ParamType*>(aligned_alloc(AlignmentLength, sizeof(ParamType)*fftSize));
        secondTransformImag = reinterpret_cast<ParamType*>(aligned_alloc(AlignmentLength, sizeof(ParamType)*fftSize));
        calculateRootsOfUnity();
    }
    
}

template<class ParamType>
void AutoCov<ParamType>::calculateRootsOfUnity()
{
    int rootsSize = (fftSize/2);
    for(int i=0; i<rootsSize; ++i)
    {
        ParamType angle = TwoPi*static_cast<ParamType>(i)/static_cast<ParamType>(fftSize);
        ParamType sinPart = std::sin(angle);
        rootsOfUnityReal[i] = std::cos(angle);
        posRootsOfUnityImag[i] = sinPart;
        negRootsOfUnityImag[i] = -sinPart;
    }
}

template<class ParamType>
int AutoCov<ParamType>::findNextPowerOfTwo(int chainLength) const
{
    //first check how many points we are using
    return static_cast<int>(std::ceil(std::log2(chainLength)));
}


}
}
}
#endif  //MCMC_ANALYSIS_AUTOCORRCALC_H
