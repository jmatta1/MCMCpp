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
#ifndef MCMC_THREADING_THREADCONTROLLER_H
#define MCMC_THREADING_THREADCONTROLLER_H
// includes for C system headers
// includes for C++ system headers
#include<atomic>
#include<thread>
#include<mutex>
#include<condition_variable>
#include<chrono>
// includes from other libraries
// includes from MCMC


namespace MarkovChainMonteCarlo
{
namespace Threading
{

/*!
 * \brief The WorkerStatus enum is used to set the state of the worker threads
 */
enum class WorkerStatus {ProcessBlack, ///<If a worker thread wakes and sees this status it updates the first half of the walkers using the second half
                         ProcessRed,  ///<If a worker thread wakes and sees this status it updates the second half of the walkers using the first half
                         Wait, ///<If a worker thread wakes and sees this status then it waits on the condition variable unti its next wait
                         Terminate ///<If a worker thread wakes and sees this status it takes steps to cease its own execution
                        };

enum class MainStatus {Wait, ///<If the main thread sees this status, it waits to be woken by the last worker at the end of a step
                       StepEnd ///<If the main thread sees this, it performs the end of step actions, and either starts the next step or stops the MCMC
                      };

/**
 * @class ThreadController
 * @ingroup Threading
 * @brief The controller object for the parallel ensemble's threads
 * @author James Till Matta
 * 
 */
class ThreadController
{
public:
    ThreadController(int numThread);
    
    void mainWaitForNewState();
    void setMainToWait();
    void setMainToStepEnd();

    bool workerSetWorkDone(){return ((++finishedWorkers) == threadCount);}
    void workerWaitForNewState();
    void setSetWorkersToProcessBlock();
    void setSetWorkersToProcessRed();
    void setSetWorkersToProcessWait();
    void setSetWorkersToProcessTerminate();
    
private:
    int threadCount; ///<Number of worker threads

    std::atomic_int terminateAcks; ///<The number of threads that have acknowledged the command to terminate
    std::atomic_int waitAck; ///<The number of threads that have acknowledged the command to wait
    std::atomic_int processRedAck;  ///<The number of threads that have acknowledged the command to process the red set
    std::atomic_int processBlackAck; ///<The number of threads that have acknowledged the command to process the black set
    
    std::atomic<MainStatus> mainState; ///<Stores the state of the main thread
    std::mutex mainWaitMutex; ///<Wait mutex for the main thread that starts everything
    std::condition_variable mainWaitCondition; ///<Wait condition variable for the main thread

    std::atomic<WorkerStatus> workerState; ///<Stores the state that the worker threads should be in
    std::mutex workerWaitMutex; ///<Wait mutex for worker threads
    std::condition_variable workerWaitCondition; ///<Wait condition variable for worker threads
    
    std::atomic_int finishedWorkers = 0; ///<Number of worker threads that have finished their work
};

}
}
#endif  //MCMC_THREADING_THREADCONTROLLER_H
