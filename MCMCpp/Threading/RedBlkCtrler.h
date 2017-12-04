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
enum class WorkerStatus {ProcessRedSave,  ///<Tells the threads to process the red batch one step and save the point to the chain
                         ProcessRedSkip,  ///<Tells the threads to process the red batch one step and skip writing
                         ProcessBlackSave, ///<Tells the threads to process the black batch one step and save the point to the chain
                         ProcessBlackSkip, ///<Tells the threads to process the black batch one step and skip writing
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
        chain(chainRef), stepAction(ptr), threadCount(numThreads), lastThread(numThreads-1), ctrlMutex(),
        ctrlWait(), midStepMutex(), midStepWait(), wrkrMutex(), wrkrWait(){}
    
    //functions for the controller thread
    /*!
     * \brief runSampling puts the thread calling this to sleep, sets up the sampling parameters, and wakes the worker threads
     * \param numSteps the number of steps to save to the chain
     * \param skipPoints the number of points to skip before saving a point to the chain
     */
    void runSampling(int numSteps, int skipPoints=0);
    
    //functions for the worker threads
    /*!
     * \brief workerWaitForTask Used by workers to wait for the worker status to change from something other than wait
     */
    void workerWaitForTask();
    
    /*!
     * \brief workerMidStepWait Used by workers to wait for the worker status to change from ProcessRedSave or ProcessRedSkip
     */
    void workerMidStepWait();
    
    /*!
     * \brief workerMidStepWait Used by workers to wait for the worker status to change from ProcessBlackSave or ProcessBlackSkip
     */
    void workerEndStepWait();
    
private:
    Chain::Chain<ParamType>& chain; ///<The chain reference for incrementing and getting itterators for the end of step action
    EndOfStepAction* stepAction = nullptr; ///<Pointer to the end of step action object
    int threadCount; ///<Number of worker threads
    int lastThread; ///<Number of threads minus one, for convenience
    
    int stepsTaken = 0; ///<Number of saved steps taken in this sampling run
    int stepsToTake = 0; ///<Number of saved steps to be taken in this sampling run, set at start of run by control thread
    int stepsSkipped = 0; ///<Number of steps skipped in this round of the sampling run
    int stepsToSkip = 0; ///<Number of steps to skip per round of the sampling run, set at start of run by control thread

    //Main thread control
    MainStatus ctrlStatus = MainStatus::Continue; ///<Status for the controller thread
    std::mutex ctrlMutex; ///<Mutex for the controller thread
    std::condition_variable ctrlWait; ///<Wait condition for the controller thread
    
    //worker control
    WorkerStatus wrkrStatus = WorkerStatus::Wait; ///<Status/command for the worker threads
    int numWorkersAtMainWait = 0; ///<Number of workers at the main wait location
    std::mutex wrkrMutex; ///<Mutex for workers to wait at the beginning / end of a step
    std::condition_variable wrkrWait; ///<Wait condition variable to hold workers at the beginning / end of a step
    
    //worker mid step control
    int numWorkersAtMidStep = 0; ///<Number of workers that are waiting at the mid step point
    std::mutex midStepMutex; ///<Mutex for workers to wait at the mid step point
    std::condition_variable midStepWait; ///<condition variable for workers to wait at the mid step point
    
    //worker end step control
    int numWorkersAtEndStep = 0; ///<Number of workers that are waiting at the mid step point
    std::mutex endStepMutex; ///<Mutex for workers to wait at the mid step point
    std::condition_variable endStepWait; ///<condition variable for workers to wait at the mid step point
};

template<class ParamType, class EndOfStepAction>
void RedBlackCtrler<ParamType, EndOfStepAction>::runSampling(int numSteps, int skipPoints)
{
    if(numSteps != 0)
    {
        //lock the control mutex
        std::unique_lock<std::mutex> ctrlLock(ctrlMutex);
        //change the control status and set up the sampling parameters
        ctrlStatus = MainStatus::Wait;
        stepsTaken = 0;
        stepsToTake = numSteps;
        stepsSkipped = 0;
        stepsToSkip = skipPoints;
        //
        if(stepsToSkip == 0)
        {
            wrkrStatus = WorkerStatus::ProcessRedSave;
        }
        else
        {
            wrkrStatus = WorkerStatus::ProcessRedSkip;
        }
        midStepWait.notify_all();
        
        while(ctrlStatus == MainStatus::Wait)
        {
            ctrlWait.wait(ctrlLock);
        }
    }
}

template<class ParamType, class EndOfStepAction>
void RedBlackCtrler<ParamType, EndOfStepAction>::workerWaitForTask()
{
    std::unique_lock<std::mutex> lock(wrkrMutex);
    while(wrkrStatus == WorkerStatus::Wait)
    {
        wrkrWait.wait(lock);
    }
}

}
}
#endif  //MCMC_THREADING_REDBLKCTRLER_H
