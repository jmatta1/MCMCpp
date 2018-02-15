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
#ifndef MCMCPP_ANALYSIS_DETAIL_AUTOCOV_H
#define MCMCPP_ANALYSIS_DETAIL_AUTOCOV_H
// includes for C system headers
// includes for C++ system headers
#include<algorithm>//for transform
#include<complex>//needed for the FFT and iFFT
#include<cmath>//needed for ceiling and log 2
// includes from other libraries
// includes from MCMCpp

namespace MCMC
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
     * \brief Deleted copy constructor
     */
    AutoCov(const AutoCov<ParamType>& rhs) = delete;
    
    /*!
     * \brief Deleted assignment operator
     */
    AutoCov<ParamType>& operator=(const AutoCov<ParamType>& rhs) = delete;
    
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
    std::complex<ParamType>* negRootsOfUnity = nullptr; ///<Stores e^(-i*2*n/N) twiddle factors used in fft
    std::complex<ParamType>* posRootsOfUnity = nullptr; ///<Stores e^(i*2*pi*k*n/N) twiddle factors used in ifft
    
    std::complex<ParamType>* firstTransform = nullptr; ///<Location where the chain is initially copied to undergo its first transform (and for its conversion to a power spectrum)
    std::complex<ParamType>* secondTransform = nullptr; ///<Location where the chain undergoes its second transform
    
    const ParamType TwoPi = (static_cast<ParamType>(2)*(static_cast<ParamType>(PiNumerator)/static_cast<ParamType>(PiDenomenator)));///<Storage of 2*pi
};

template<class ParamType>
AutoCov<ParamType>::~AutoCov()
{
    if(negRootsOfUnity != nullptr) delete[] negRootsOfUnity;
    if(posRootsOfUnity != nullptr) delete[] posRootsOfUnity;
    if(firstTransform != nullptr) delete[] firstTransform;
    if(secondTransform != nullptr) delete[] secondTransform;
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
    ParamType normValue = secondTransform[0].real();
    std::transform(secondTransform, secondTransform+chainLength, output,
                   [&normValue](const std::complex<ParamType>& v) -> ParamType {return (v.real()/normValue);});
}

template<class ParamType>
void AutoCov<ParamType>::performInverseFft()
{
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
                std::complex<ParamType> temp1 = (posRootsOfUnity[twiddleIndex] * secondTransform[k+j+m2]);
                std::complex<ParamType> temp2 = secondTransform[k+j];
                secondTransform[k+j] = (temp2 + temp1);
                secondTransform[k+j+m2] = (temp2 - temp1);
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
        std::complex<ParamType> temp = firstTransform[i];
        secondTransform[bitReverse(i)] = (temp.real()*temp.real()+temp.imag()*temp.imag());
    }
}

template<class ParamType>
void AutoCov<ParamType>::performForwardFft()
{//simple implementation of Cooley-Tukey
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
                std::complex<ParamType> temp1 = (negRootsOfUnity[twiddleIndex] * firstTransform[k+j+m2]);
                std::complex<ParamType> temp2 = firstTransform[k+j];
                firstTransform[k+j] = (temp2 + temp1);
                firstTransform[k+j+m2] = (temp2 - temp1);
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
        firstTransform[bitReverse(i)] = (chain[i] - avg);
    }
    
    for(int i=chainLength; i<fftSize; ++i)
    {
        firstTransform[bitReverse(i)] = static_cast<ParamType>(0);
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
        if(negRootsOfUnity != nullptr) delete[] negRootsOfUnity;
        if(posRootsOfUnity != nullptr) delete[] posRootsOfUnity;
        if(firstTransform  != nullptr) delete[] firstTransform;
        if(secondTransform != nullptr) delete[] secondTransform;
        negRootsOfUnity = new std::complex<ParamType>[fftSize/2];
        posRootsOfUnity = new std::complex<ParamType>[fftSize/2];
        firstTransform = new std::complex<ParamType>[fftSize];
        secondTransform = new std::complex<ParamType>[fftSize];
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
        ParamType cosPart = std::cos(angle);
        ParamType sinPart = std::sin(angle);
        negRootsOfUnity[i] = std::complex<ParamType>{cosPart, -sinPart};
        posRootsOfUnity[i] = std::complex<ParamType>{cosPart, sinPart};
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
#endif  //MCMCPP_ANALYSIS_AUTOCORRCALC_H
