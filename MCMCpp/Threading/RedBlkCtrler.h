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


namespace MCMC
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
    RedBlackCtrler(int numThreads, Chain::Chain<ParamType>& chainRef, EndOfStepAction* ptr = nullptr):
        threadCount(numThreads), lastThread(numThreads-1), chain(chainRef), stepAction(ptr){}
    
    //functions for the controller thread
    void controllerWaitForContinue();
    void runSampling(int numSteps, int skipPoints=0);
    
    //functions for the worker threads
    void workerWaitForChange();
    void workerMidStepWait();
    void workerEndStepWait();
    
private:
    int threadCount; ///<Number of worker threads
    Chain::Chain<ParamType>& chain; ///<The chain reference for incrementing and getting itterators for the end of step action
    EndOfStepAction* stepAction = nullptr; ///<Pointer to the end of step action object
    
    int stepsTaken = 0; ///<Number of saved steps taken in this sampling run
    int stepsToTake; ///<Number of saved steps to be taken in this sampling run
    int stepsSkipped = 0; ///<Number of steps skipped in this round of the sampling run
    int stepsToSkip; ///<Number of steps to skip per round of the sampling run

    MainStatus ctrlStatus; ///<Status for the controller thread
    std::mutex ctrlMutex; ///<Mutex for the controller thread
    std::condition_variable ctrlWait; ///<Wait condition for the controller thread
    
    WorkerStatus wrkrStatus; ///<Status/command for the worker threads
    int numWorkersAtMidStep = 0; ///<Number of workers that are waiting at the mid step point
    std::mutex midStepMutex; ///<Mutex for workers to wait at the mid step point
    std::condition_variable midStepWait; ///<condition variable for workers to wait at the mid step point
    int numWorkersAtMainWait = 0; ///<Number of workers at the main wait location
    std::mutex wrkrMutex; ///<Mutex for workers to wait at the beginning / end of a step
    std::condition_variable wrkrWait; ///<Wait condition variable to hold workers at the beginning / end of a step
};

}
}
#endif  //MCMC_THREADING_REDBLKCTRLER_H
