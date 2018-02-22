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
#ifndef MCMCPP_THREADING_REDBLKCTRLER_H
#define MCMCPP_THREADING_REDBLKCTRLER_H
// includes for C system headers
// includes for C++ system headers
#include<mutex>
#include<condition_variable>
#include<atomic>
#include<cassert>
#include<thread>
#include<chrono>
// includes from other libraries
// includes from MCMCpp
#include"../Chain/Chain.h"


namespace MCMC
{
namespace Threading
{

/*!
 * \brief The WorkerStatus enum is used to set the state of the worker threads
 */
enum class WorkerStatus {ProcessRed,  ///<Tells the threads to process the red batch forward one step
                         ProcessBlack, ///<Tells the threads to process the black batch forward one step
                         Wait, ///<Worker thread wait to be woken state
                         Terminate ///< Worker thread cease execution state
                        };

enum class MainStatus {Wait, ///<If the main thread sees this status, it waits to be woken
                       Continue ///<If the main thread sees this, it resumes its flow
                      };

/**
 * @class RedBlackCtrler
 * @ingroup Threading
 * @brief The controller object for the parallel ensemble's threads
 * @author James Till Matta
 * 
 */
template<class ParamType, class EndOfStepAction, bool UseSpinLocks>
class RedBlackCtrler
{};

}
}
//pull in the specializations
#include"RedBlkCtrlerNormalLock.h"
#include"RedBlkCtrlerSpinLock.h"

#endif  //MCMCPP_THREADING_REDBLKCTRLER_H
