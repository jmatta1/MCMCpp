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
#ifndef MCMC_UTILITY_AUTOCORRCALC_H
#define MCMC_UTILITY_AUTOCORRCALC_H

namespace MarkovChainMonteCarlo
{
namespace Utility
{
/*!
 * @class AutoCorrCalc
 * @ingroup Utility
 * @brief A class to calculate the autocorrelation times of the parameters that were MCMCed
 * @author James Till Matta
 * 
 * @tparam ParamType The floating point type for which this calculation will be carried out
 * @tparam BlockSize Number of steps to store in a block of the linked list that stores the chain 
 */
template<class ParamType, int BlockSize>
class AutoCorrCalc
{
    /*!
     * \brief AutoCorrCalc constructs a new AutoCorrCalc object
     * \param numParams The number of parameters in each sample
     * \param numWalkers The number of walkers in the ensemble
     */
    AutoCorrCalc(int numParams, int numWalkers);
    /*!
     * \brief setAutoCorrelationTimeParameters
     * \param minAutoCorrTimes The minimum required number of autocorrelation times the algorithm needs to examine
     * \param step The increase in window size for each iteration of the algorithm
     * \param loWin The minimum window size
     * \param hiWin The maximum window size
     * \param fast If true, only use the first power of two samples to accellerate calculation using FFT
     */
    void setAutoCorrelationTimeParameters(int minAutoCorrTimes=10, int step=1, int loWin=10, int hiWin=10000, bool fast=false);
    
    /*!
     * \brief getAutoCorrelationTime Calculates the correlation time for a given parameter
     * \param paramIndex The index of the parameter [0, numParameter)
     * \return The autocorrelation time in samples for parameter # paramIndex
     */
    ParamType getAutoCorrelationTime(int paramIndex);
    
private:
    
};

}
}
#endif  //MCMC_UTILITY_AUTOCORRCALC_H
