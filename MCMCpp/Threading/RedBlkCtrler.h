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
#ifndef MCMC_THREADING_REDBLKCTRLER_H
#define MCMC_THREADING_REDBLKCTRLER_H
// includes for C system headers
// includes for C++ system headers
#include<atomic>
#include<thread>
#include<mutex>
#include<condition_variable>
#include<chrono>
// includes from other libraries
// includes from MCMC
#include"../Chain/Chain.h"


namespace MarkovChainMonteCarlo
{
namespace Threading
{

/*!
 * \brief The WorkerStatus enum is used to set the state of the worker threads
 */
enum class WorkerStatus {ProcessBlack, ///<
                         ProcessRed,  ///<
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
template<class ParamType, class EndOfStepAction>
class RedBlackCtrler
{
public:
    RedBlackCtrler(int numThreads, Chain::Chain<ParamType>& chainRef,
                   EndOfStepAction* ptr = nullptr):threadCount(numThreads){}
    
    //functions for the controller thread
    void controllerWaitForContinue();
    void runSampling(int numSteps, int skipPoints=0);
    
    //functions for the worker threads
    void workerWaitForChange();
    void workerMidStepWait();
    void workerEndStepWait();
    
private:
    int threadCount; ///<Number of worker threads
};

}
}
#endif  //MCMC_THREADING_REDBLKCTRLER_H
