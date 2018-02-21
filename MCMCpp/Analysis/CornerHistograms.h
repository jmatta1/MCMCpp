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
#ifndef MCMCPP_ANALYSIS_CORNERHISTOGRAMS_H
#define MCMCPP_ANALYSIS_CORNERHISTOGRAMS_H
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
 * @class CornerHistograms
 * @ingroup Analysis
 * @brief A class to calculate the single and bi-axial histograms of the points in chains
 * @author James Till Matta
 * 
 * @tparam ParamType The floating point type for which this calculation will be carried out
 * 
 * The user can skip parts of the chain to use every nth step in or simply use every point
 * Care should be taken in choice of binning a good rule of thumb is that the memory size
 * used by CornerHistograms is going to be approximately
 * Size in bytes = (2*n*(n-1)*b*b + 4*n*b) where n is the number of parameters
 * and b is the number of bins
 */
template<class ParamType>
class CornerHistograms
{
public:
    typedef Chain::ChainStepIterator<ParamType> IttType;
    
    /*!
     * \brief Constructor for Corner Histogram calculator
     * \param numParams The number of parameters per point
     * \param numWalkers The number of walkers per step in the chain
     * \param binsPerAxis The desired number of bins per axis in the histograms, 100 bins per axis is a somewhat sane default value
     */
    CornerHistograms(int numParams, int numWalkers, int binsPerAxis=100):
        pCount(numParams), wCount(numWalkers), bCount(binsPerAxis)
    {
        assert(pCount > 0);
        assert(wCount > 0);
        assert(bCount > 1);
        numTwoAxis = ((pCount*(pCount-1))/2);
        //allocate the array of pointers two the 2d hists
        twoAxisHists = new int*[numTwoAxis];
        //allocate each of the 2D hists
        for(int i=0; i<numTwoAxis; ++i)
        {
            twoAxisHists[i] = new int[bCount*bCount];
        }
        //allocate the array of 1D hists
        singleAxisHists = new int[bCount*pCount];
        //allocate the set of parameter bounds and set them to extreme values
        paramBounds = new ParamType[2*pCount];
    }
    
    ~CornerHistograms()
    {
        delete[] singleAxisHists;
        for(int i=0; i<numTwoAxis; ++i)
        {
            delete[] twoAxisHists[i];
        }
        delete[] twoAxisHists;
    }

    /*!
     * \brief Histograms the points in the chain into all concievable 1 and 2D hisograms
     * \param start The iterator at the start of the chain
     * \param end The iterator at the end of the chain
     * \param sliceInterval The relative index of sample to use, must be >= 1. a value of 1 uses every sample, 2 uses every other sample, so on and so forth
     */
    void calculateHistograms(IttType start, IttType end, int sliceInterval=1)
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
    }
    
    /*!
     * \brief Saves all the histograms in a standard csv format
     * \param fileNameBase The root of the file name to which specifiers are written
     */
    void saveHistsCsvFormat(const std::string& fileNameBase)
    {
        //iterate through the parameter pair combinations
        for(int i=0; i<pCount; ++i)
        {//first write the 1D histogram
            write1dCsv(i, fileNameBase);
            //then handle the dual parameter histograms
            int cornerOffset = ((i*(i-1))/2);
            for(int j=0; j<i; ++j)
            {
                write2dCsv(i, j, fileNameBase);
            }
        }
    }
    
    /*!
     * \brief Returns the low edge for the bin number specified for a given parameter number
     * \param pNum The parameter number whose binning is to be examined
     * \param binNum The bin number to get the low edge for
     * \return The low edge of bin number binNum for parameter pNum
     */
    ParamType getHistBinLowEdge(int pNum, int binNum)
    {
        return (paramBounds[2*pNum] + (binNum*paramBounds[2*pNum+1]));
    }
    
    /*!
     * \brief Returns the low edge for the bin number specified for a given parameter number
     * \param pNum The parameter number whose binning is to be examined
     * \param binNum The bin number to get the low edge for
     * \return The low edge of bin number binNum for parameter pNum
     */
    ParamType getHistBinHighEdge(int pNum, int binNum)
    {
        return (paramBounds[2*pNum] + ((binNum+1)*paramBounds[2*pNum+1]));
    }
    
    /*!
     * \brief Retrieve the specified bin for parameter number p1
     * \param pNum The parameter whose bin is to be examined
     * \param bin The bin number of the histogram to get
     * \return The value of the histogram for p1 for bin
     */
    int get1dHistBin(int pNum, int bin)
    {
        return singleAxisHists[pNum*bCount+bin];
    }
    
    /*!
     * \brief Retrieve the specified 2D histogram bin
     * \param p1 The first parameter whose 2D histogram is to be examined must be greater than p2
     * \param p2 The second parameter whose 2D histogram is to be examined must be less than p1
     * \param binx The xaxis bin number to retrieve from the histogram 
     * \param biny The yaxis bin number to retrieve from the histogram
     * \return The value of the 2d histogram for p1 and p2 for binx and biny
     */
    void get2dHistBin(int p1, int p2, int binx, int biny)
    {
        int cornerOffset = ((p1*(p1-1))/2);
        return twoAxisHists[cornerOffset + p2][biny*bCount+binx];
    }
    
private:
    /*!
     * \brief Write the specified 2D parameter histogram into csv format
     * \param p1 the parameter of the X-axis
     * \param p2 the parameter of the Y-axis
     * \param fileNameBase The base filename from which the output filename will be calculated
     */
    void write2dCsv(int p1, int p2, const std::string& fileNameBase)
    {
        std::ostringstream namer;
        namer<<fileNameBase<<"_p"<<p1<<"_p"<<p2<<".csv";
        std::ofstream out;
        out.open(namer.str().c_str());
        out<<"# Lines starting with a '#' in the first column are ignored\n";
        out<<"# X-axis: nbins, first bin low edge, last bin high edge\n";
        //write the x-axis parameters
        out<<bCount<<", "<<paramBounds[2*p1]<<", "<<(paramBounds[2*p1] + (bCount*paramBounds[2*p1+1]))<<"\n";
        out<<"# Y-axis: nbins, first bin low edge, last bin high edge\n";
        //write the y-axis parameters
        out<<bCount<<", "<<paramBounds[2*p2]<<", "<<(paramBounds[2*p2] + (bCount*paramBounds[2*p2+1]))<<"\n";
        out<<"# x-bin number, y-bin number, value\n";
        int cornerOffset = ((p1*(p1-1))/2);
        int* hist = twoAxisHists[cornerOffset+p2];
        for(int i=0; i<bCount; ++i)
        {
            int offset = i*bCount;
            for(int j=0; j<bCount; ++j)
            {
                out<<i<<", "<<j<<", "<<hist[offset+j]<<std::endl;
            }
        }
        out<<std::flush;
        out.close();
    }
    
    /*!
     * \brief Writes the specified single parameter histogram into csv format
     * \param histNum The parameter number of the 1D histogram
     * \param fileNameBase The base filename from which the output filename will be calculated
     */
    void write1dCsv(int histNum, const std::string& fileNameBase)
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
            out<<i<<", "<<singleAxisHists[offset+i]<<"\n";
        }
        out<<std::flush;
        out.close();
    }
    
    /*!
     * \brief Takes a point and increments the histograms as needed using the point
     * \param pt The point to be used to increment histograms
     */
    void handlePoint(ParamType* pt)
    {
        //iterate through the parameter pair combinations
        for(int i=0; i<pCount; ++i)
        {
            //first handle the single parameter histogram
            int iBin = ((pt[i] - paramBounds[2*i])/paramBounds[2*i+1]);
            ++(singleAxisHists[i*bCount+iBin]);
            //then handle the dual parameter histograms
            int cornerOffset = ((i*(i-1))/2);
            for(int j=0; j<i; ++j)
            {
                int jBin = ((pt[j]-paramBounds[2*j])/paramBounds[2*j+1]);
                ++(twoAxisHists[cornerOffset+j][jBin*bCount+iBin]);
            }
        }
    }
    
    /*!
     * \brief Finds the minimum and maximum value of every parameter
     * \param start The step iterator at the start of the chain
     * \param end The step iterator at the end of the chain
     * \param sliceInterval The relative index of sample to use, must be >= 1. a value of 1 uses every sample, 2 uses every other sample, so on and so forth
     */
    void findBinning(IttType start, IttType end, int sliceInterval)
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
    
    /*!
     * \brief Checks if each parameter constitutes an extreme (min or max) for a given point
     * \param pt The point whose values are to be tested
     */
    void testPtExtremes(ParamType* pt)
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

    /*!
     * \brief Resets all the storage in the class to zero, or relevant value for starting
     */
    void zeroStorage()
    {
        int histSize = (bCount*bCount);
        //zero each of the 2D hists
        for(int i=0; i<numTwoAxis; ++i)
        {
            for(int j=0; j<histSize; ++j)
            {
                twoAxisHists[i][j] = 0;
            }
        }
        //zero all the 1D hists
        int arrSize = bCount*pCount;
        for(int i=0; i<arrSize; ++i)
        {
            singleAxisHists[i] = 0;
        }
        //set the parameter bounds to their initial positions
        for(int i=0; i<pCount; ++i)
        {
            paramBounds[2*i] = std::numeric_limits<ParamType>::max();
            paramBounds[2*i+1] = std::numeric_limits<ParamType>::min();
        }
    }
    
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
    
    int** twoAxisHists = nullptr; ///<Holds the list of pointers to two-d histograms
    int* singleAxisHists = nullptr; ///<Holds the set of 1-d histograms, one per row of the 2d array
    ParamType* paramBounds = nullptr; ///<Holds the min and max of every parameter
    int pCount;///<Holds the number of parameters per point
    int wCount;///<Holds the number of walkers per step
    int bCount;///<Holds the number of bins per axis
    int numTwoAxis = 0;///<Holds the number of two-axis histograms
};

}
}

#endif  //MCMCPP_ANALYSIS_CORNERHISTOGRAMS_H
