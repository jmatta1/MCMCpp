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
#include<mutex>
#include<condition_variable>
#include<atomic>
#include<cassert>
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
template<class ParamType, class EndOfStepAction>
class RedBlackCtrler
{
public:
    RedBlackCtrler(int numThreads, Chain::Chain<ParamType>& chainRef, EndOfStepAction* ptr = nullptr):
        chain(chainRef), stepAction(ptr), threadCount(numThreads), ctrlMutex(),
        ctrlWait(), midStepWait(), wrkrMutex(), wrkrWait(){}
    
    //functions for the controller thread
    /*!
     * \brief runSampling puts the thread calling this to sleep, sets up the sampling parameters, and wakes the worker threads
     * \param numSteps the number of steps to save to the chain
     * \param skipPoints the number of points to skip before saving a point to the chain
     */
    void runSampling(int numSteps, int skipPoints=0);
    
    /*!
     * \brief getNumWorkersWaiting returns the number of worker threads waiting for work
     * \return The number of threads waiting on the main condition variable
     */
    int getNumWorkersWaiting(){return numWorkersAtMainWait.load();}
    
    /*!
     * \brief terminateWorkers Sets the status of all workers to terminate and ensures that all worker condition variables are notified
     */
    void terminateWorkers();
    
    /*!
     * \brief getTerminatedWorkerCount Returns the number of worker threads that have terminated
     * \return The number of threads that have terminated
     */
    int getTerminatedWorkerCount(){return wrkrTerm;}
    
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
    
    /*!
     * \brief workerAckTerminate Used by the worker threads to signal that they have exited their event loop
     */
    void workerAckTerminate(){wrkrTerm.fetch_add(1);}
    
    /*!
     * \brief workerGetStatus Returns the current state the worker should be in
     * \return The currently set worker state
     */
    WorkerStatus workerGetStatus(){return wrkrStatus.load();}
    /*!
     * \brief workerGetSavingPoint Returns if the worker should be saving points or not
     * \return The current setting for saving points
     */
    bool workerGetSavingPoint(){return savePoints.load();}
    
private:
    Chain::Chain<ParamType>& chain; ///<The chain reference for incrementing and getting itterators for the end of step action
    EndOfStepAction* stepAction = nullptr; ///<Pointer to the end of step action object
    int threadCount; ///<Number of worker threads
    
    int stepsTaken = 0; ///<Number of saved steps taken in this sampling run
    int stepsToTake = 0; ///<Number of saved steps to be taken in this sampling run, set at start of run by control thread
    int stepsSkipped = 0; ///<Number of steps skipped in this round of the sampling run
    int stepsToSkip = 0; ///<Number of steps to skip per round of the sampling run, set at start of run by control thread

    //Main thread control
    std::atomic<MainStatus> ctrlStatus = MainStatus::Continue; ///<Status for the controller thread
    std::mutex ctrlMutex; ///<Mutex for the controller thread
    std::condition_variable ctrlWait; ///<Wait condition for the controller thread
    
    //worker control
    std::atomic<WorkerStatus> wrkrStatus = WorkerStatus::Wait; ///<Status/command for the worker threads
    std::atomic_bool savePoints = true; ///<Whether or not worker threads should tell their walkers to save this step
    std::atomic_int numWorkersAtMainWait = 0; ///<Number of workers at the main wait location, queryable from the outside, so it is atomic
    std::mutex wrkrMutex; ///<Mutex for workers to wait at the beginning / end of a step
    std::condition_variable wrkrWait; ///<Wait condition variable to hold workers when they are not stepping
    std::atomic_int wrkrTerm = 0;
    
    //worker mid step control
    int numWorkersAtMidStep = 0; ///<Number of workers that are waiting at the mid step point
    std::condition_variable midStepWait; ///<condition variable for workers to wait at the mid step point
    
    //worker end step control
    int numWorkersAtEndStep = 0; ///<Number of workers that are waiting at the mid step point
    std::condition_variable endStepWait; ///<condition variable for workers to wait at the mid step point
};

template<class ParamType, class EndOfStepAction>
void RedBlackCtrler<ParamType, EndOfStepAction>::terminateWorkers()
{
    std::unique_lock<std::mutex> ctrlLock(ctrlMutex);
    std::unique_lock<std::mutex> wrkrLock(wrkrMutex);
    wrkrTerm.store(0);
    wrkrStatus.store(WorkerStatus::Terminate);
    wrkrWait.notify_all();
    midStepWait.notify_all();
    endStepWait.notify_all();
}

template<class ParamType, class EndOfStepAction>
void RedBlackCtrler<ParamType, EndOfStepAction>::runSampling(int numSteps, int skipPoints)
{
    //Check for critical errors
    //Where Either the number of steps to sample or the number of points to skip was less than 1
    assert(numSteps <= 0);
    assert(skipPoints <= 0);
    //lock the control thread mutex
    std::unique_lock<std::mutex> ctrlLock(ctrlMutex);
    
    {//artificial scope to force destruction of lock on worker mutex when done modifying variables that workers can use
        //lock the worker mutex to prevent unforseen races
        std::unique_lock<std::mutex> wrkrLock(wrkrMutex);
        //change the control status and set up the sampling parameters
        ctrlStatus.store(MainStatus::Wait);
        stepsTaken = 0;
        stepsToTake = numSteps;
        stepsSkipped = 0;
        stepsToSkip = skipPoints;
        //if we are not doing sample slicing, jump straight to a save step
        //otherwise set to a skip step
        if(stepsToSkip == 0)
        {
            savePoints.store(true);
        }
        else
        {
            savePoints.store(false);
        }
        //Set the number of threads at the midstep checkpoint to zero
        numWorkersAtMidStep = 0;
        //Set the processing mode to start stepping the ensemble
        wrkrStatus.store(WorkerStatus::ProcessRed);
        //wake all the worker threads sitting on the primary worker condition variable
        wrkrWait.notify_all();
        //the worker mutex will be unlocked now so the worker threads can wake up
    }
    //make the control thread (this thread) sleep until it is signalled that sampling is finished
    while(ctrlStatus.load() == MainStatus::Wait)
    {
        ctrlWait.wait(ctrlLock);
    }
}

template<class ParamType, class EndOfStepAction>
void RedBlackCtrler<ParamType, EndOfStepAction>::workerWaitForTask()
{
    //lock the worker mutex, increment the wait counter, then wait
    std::unique_lock<std::mutex> lock(wrkrMutex);
    numWorkersAtMainWait.fetch_add(1);
    while(wrkrStatus.load() == WorkerStatus::Wait)
    {
        wrkrWait.wait(lock);
    }
}

template<class ParamType, class EndOfStepAction>
void RedBlackCtrler<ParamType, EndOfStepAction>::workerMidStepWait()
{
    //Note: if we are at the mid step the status is ProcessRed, always
    std::unique_lock<std::mutex> lock(wrkrMutex);
    ++numWorkersAtMidStep;
    if(numWorkersAtMidStep == threadCount)
    {//we are the last thread to finish:
        //reset the number threads waiting at the next checkpoint, change worker state, and wake the workers
        numWorkersAtEndStep = 0;
        wrkrStatus.store(WorkerStatus::ProcessBlack);
        midStepWait.notify_all();
    }
    else
    {//we are not the last thread to finish, so wait until the others are ready
        while(wrkrStatus.load() == WorkerStatus::ProcessRed)
        {
            midStepWait.wait(lock);
        }
    }
    
}

template<class ParamType, class EndOfStepAction>
void RedBlackCtrler<ParamType, EndOfStepAction>::workerEndStepWait()
{
    //Note: if we are at the end step the status is ProcessBlack, always
    std::unique_lock<std::mutex> lock(wrkrMutex);
    ++numWorkersAtEndStep;
    if(numWorkersAtEndStep == threadCount)
    {//we are the last thread to finish the step
        //create a variable to store if we are at the end of sampling for when the time comes to update worker thread states
        bool endOfSampling = false;
        //lock the control mutex to prevent unforseen races, this should be basicly contention free and thus fast
        std::unique_lock<std::mutex> ctrlLock(ctrlMutex);
        //first increment the chain
        chain.incrementChainStep();
        //now call the post step action if there is one
        if(stepAction != nullptr)
        {
             stepAction->performAction(chain.getStepIteratorBegin(), chain.getStepIteratorEnd());
        }
        //next increment the appropriate step counter and adjust the write mode as needed
        if(savePoints.load())
        {
            ++stepsTaken;
            //check if we are done sampling (since we can only finish sampling on a written point
            if(stepsTaken == stepsToTake)
            {
                endOfSampling = true;
            }
            //check if we need to shift to skip mode for the next X points
            if(stepsToSkip != 0)
            {//whenever we write a point we need to skip the next batch before writing again
                savePoints.store(false);
                stepsSkipped = 0;
            }
        }
        else
        {
            ++stepsSkipped;
            //we cannot be done sampling if we are in a stretch of skipping
            //check if we need to save the next step
            if(stepsToSkip == stepsSkipped)
            {
                savePoints.store(false);
            }
        }
        //handle if we are at the end of sampling
        if(endOfSampling)
        {
            numWorkersAtMainWait.store(0);
            wrkrStatus.store(WorkerStatus::Wait);
            ctrlStatus.store(MainStatus::Continue);
            endStepWait.notify_all();
            ctrlWait.notify_all();
        }
        else
        {
            numWorkersAtMidStep = 0;
            wrkrStatus.store(WorkerStatus::ProcessRed);
            endStepWait.notify_all();
        }
    }
    else
    {//we are not the last thread, wait for further instructions
        while(wrkrStatus.load() == WorkerStatus::ProcessBlack)
        {
            endStepWait.wait(lock);
        }
    }
}

}
}
#endif  //MCMC_THREADING_REDBLKCTRLER_H
