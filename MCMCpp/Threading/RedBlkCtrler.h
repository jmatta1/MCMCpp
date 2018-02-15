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
template<class ParamType, class EndOfStepAction>
class RedBlackCtrler
{
public:
    /*!
     * \brief RedBlackCtrler Constructor for the controller class that acts as the central switchyard for the worker threads
     * \param numThreads the number of threads to be controlled
     * \param chainRef A reference to the chain object
     * \param ptr A pointer to the EndOfStepAction class, nullptr if there will be no end of step action
     */
    RedBlackCtrler(int numThreads, Chain::Chain<ParamType>& chainRef, EndOfStepAction* ptr = nullptr):
        chain(chainRef), stepAction(ptr), threadCount(numThreads),
        lastThread(numThreads-1), stepsTaken(0), stepsSkipped(0),
        ctrlStatus(MainStatus::Continue), ctrlMutex(), ctrlWait(),
        wrkrStatus(WorkerStatus::Wait), savePoints(true), wrkrIndex(0),
        numWorkersAtMainWait(0), wrkrMutex(), wrkrWait(), wrkrTerm(0),
        numWorkersAtMidStep(0), numWorkersAtEndStep(0){}
    
    //functions for the controller thread
    /*!
     * \brief runSampling puts the thread calling this to sleep, sets up the sampling parameters, and wakes the worker threads
     * \param numSteps the number of steps to save to the chain
     * \param skipPoints the number of points to skip before saving a point to the chain
     */
    void runSampling(int numSteps, int skipPoints=0);
    
    /*!
     * \brief samplingComplete Returns true if the requested number of samples has been taken, otherwise false
     * \return True if the number of steps taken is equal to the number of steps requested, false otherwise
     */
    bool samplingComplete(){return (stepsTaken.load(std::memory_order_acquire) == stepsToTake);}
    
    /*!
     * \brief getNumWorkersWaiting returns the number of worker threads waiting for work
     * \return The number of threads waiting on the main condition variable
     */
    int getNumWorkersWaiting(){return numWorkersAtMainWait.load(std::memory_order_acquire);}
    
    /*!
     * \brief terminateWorkers Sets the status of all workers to terminate and ensures that all worker condition variables are notified
     */
    void terminateWorkers();
    
    /*!
     * \brief getTerminatedWorkerCount Returns the number of worker threads that have terminated
     * \return The number of threads that have terminated
     */
    int getTerminatedWorkerCount(){return wrkrTerm.load(std::memory_order_acquire);}
    
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
    void workerAckTerminate();
    
    /*!
     * \brief workerGetStatus Returns the current state the worker should be in
     * \return The currently set worker state
     */
    WorkerStatus workerGetStatus(){return wrkrStatus.load(std::memory_order_acquire);}
    
    /*!
     * \brief workerGetSavingPoint Returns if the worker should be saving points or not
     * \return The current setting for saving points
     */
    bool workerGetSavingPoint(){return savePoints.load(std::memory_order_acquire);}
    
    /*!
     * \brief workerGetWorkIndex Gets the next work index for the worker thread
     * \return the next worker index
     */
    int workerGetWorkIndex(){return wrkrIndex.fetch_add(1, std::memory_order_acq_rel);}
    
private:
    Chain::Chain<ParamType>& chain; ///<The chain reference for incrementing and getting itterators for the end of step action
    EndOfStepAction* stepAction = nullptr; ///<Pointer to the end of step action object
    int threadCount; ///<Number of worker threads
    int lastThread; ///<Number of worker threads - 1, stored for convenience
    
    std::atomic_int stepsTaken; ///<Number of saved steps taken in this sampling run
    int stepsToTake = 0; ///<Number of saved steps to be taken in this sampling run, set at start of run by control thread
    std::atomic_int stepsSkipped; ///<Number of steps skipped in this round of the sampling run
    int stepsToSkip = 0; ///<Number of steps to skip per round of the sampling run, set at start of run by control thread
    
    //Main thread control
    std::atomic<MainStatus> ctrlStatus; ///<Status for the controller thread
    std::mutex ctrlMutex; ///<Mutex for the controller thread
    std::condition_variable ctrlWait; ///<Wait condition for the controller thread
    
    //worker action control
    std::atomic<WorkerStatus> wrkrStatus; ///<Status/command for the worker threads
    std::atomic_bool savePoints; ///<Whether or not worker threads should tell their walkers to save this step
    std::atomic_int wrkrIndex; ///<The index of the next walker in the set to update
    
    //worker main wait control
    std::atomic_int numWorkersAtMainWait; ///<Number of workers at the main wait location, queryable from the outside, so it is atomic
    std::mutex wrkrMutex; ///<Mutex for workers to wait at the beginning / end of a step
    std::condition_variable wrkrWait; ///<Wait condition variable to hold workers when they are not stepping
    
    //worker terminate control
    std::atomic_int wrkrTerm; ///<Number of workers that have acknowledged the terminate command / terminated
    
    //worker mid step control
    std::atomic_int numWorkersAtMidStep; ///<Number of workers that are waiting at the mid step point
    
    //worker end step control
    std::atomic_int numWorkersAtEndStep; ///<Number of workers that are waiting at the mid step point
};

template<class ParamType, class EndOfStepAction>
void RedBlackCtrler<ParamType, EndOfStepAction>::terminateWorkers()
{
    {
        std::unique_lock<std::mutex> wrkrLock(wrkrMutex);
        wrkrTerm.store(0, std::memory_order_release);
        wrkrStatus.store(WorkerStatus::Terminate);
        wrkrWait.notify_all();
    }
    std::unique_lock<std::mutex> ctrlLock(ctrlMutex);
    while(wrkrTerm.load(std::memory_order_acquire) != threadCount)
    {
        ctrlWait.wait_for(ctrlLock, std::chrono::microseconds(10));
    }
}

template<class ParamType, class EndOfStepAction>
void RedBlackCtrler<ParamType, EndOfStepAction>::runSampling(int numSteps, int skipPoints)
{
    //Check for critical errors
    //Where Either the number of steps to sample or the number of points to skip was less than 1
    assert(numSteps > 0);
    assert(skipPoints >= 0);
    {//artificial scope to force destruction of lock on worker mutex when done modifying variables that workers can use
        //lock the worker mutex to prevent unforseen races
        std::unique_lock<std::mutex> wrkrLock(wrkrMutex);
        //change the control status and set up the sampling parameters
        ctrlStatus.store(MainStatus::Wait, std::memory_order_release);
        stepsTaken.store(0, std::memory_order_release);;
        stepsToTake = numSteps;
        stepsSkipped.store(0, std::memory_order_release);
        stepsToSkip = skipPoints;
        //if we are not doing sample slicing, jump straight to a save step
        //otherwise set to a skip step
        if(stepsToSkip == 0)
        {
            savePoints.store(true, std::memory_order_release);
        }
        else
        {
            savePoints.store(false, std::memory_order_release);
        }
        //set the starting index of threads to zero
        wrkrIndex.store(0, std::memory_order_release);
        //Set the number of threads at the midstep checkpoint to zero
        numWorkersAtMidStep.store(0, std::memory_order_release);
        //Set the processing mode to start stepping the ensemble
        wrkrStatus.store(WorkerStatus::ProcessRed, std::memory_order_release);
        //wake all the worker threads sitting on the primary worker condition variable
        wrkrWait.notify_all();
        //the worker mutex will be unlocked now so the worker threads can wake up
    }
    //lock the control thread mutex
    std::unique_lock<std::mutex> ctrlLock(ctrlMutex);
    //make the control thread (this thread) sleep until it is signalled that sampling is finished
    while(ctrlStatus.load(std::memory_order_acquire) == MainStatus::Wait)
    {
        ctrlWait.wait(ctrlLock);
    }
}

template<class ParamType, class EndOfStepAction>
void RedBlackCtrler<ParamType, EndOfStepAction>::workerWaitForTask()
{
    //lock the worker mutex, increment the wait counter, then wait
    std::unique_lock<std::mutex> lock(wrkrMutex);
    numWorkersAtMainWait.fetch_add(1, std::memory_order_acq_rel);
    while(wrkrStatus.load(std::memory_order_acquire) == WorkerStatus::Wait)
    {
        wrkrWait.wait(lock);
    }
}

template<class ParamType, class EndOfStepAction>
void RedBlackCtrler<ParamType, EndOfStepAction>::workerMidStepWait()
{
    int temp = numWorkersAtMidStep.fetch_add(1, std::memory_order_acq_rel);
    if(temp == lastThread)
    {//we are the last thread swap things around and go
        numWorkersAtEndStep.store(0, std::memory_order_release);
        wrkrIndex.store(0, std::memory_order_release);
        wrkrStatus.store(WorkerStatus::ProcessBlack, std::memory_order_release);
    }
    else
    {//otherwise, spin until the state changes to something other than ProcessRed
        while(wrkrStatus.load(std::memory_order_acquire) == WorkerStatus::ProcessRed){}
    }
    
}

template<class ParamType, class EndOfStepAction>
void RedBlackCtrler<ParamType, EndOfStepAction>::workerEndStepWait()
{
    //Note: if we are at the end step the status is ProcessBlack, always
    int temp = numWorkersAtEndStep.fetch_add(1, std::memory_order_acq_rel);
    if(temp == lastThread)
    {//we are the last thread to finish the step
        //create a variable to store if we are at the end of sampling for when the time comes to update worker thread states
        bool endOfSampling = false;
        //first increment the chain, if we have reached the max memory allowable for the chain, declare sampling at an end
        if(Chain::IncrementStatus::EndOfChain == chain.incrementChainStep())
        {
            endOfSampling = true;
        }
        //now call the post step action if there is one
        if(stepAction != nullptr)
        {
             stepAction->performAction(chain.getStepIteratorBegin(), chain.getStepIteratorEnd());
        }
        //next increment the appropriate step counter and adjust the write mode as needed
        if(savePoints.load(std::memory_order_acquire))
        {
            int temp = (stepsTaken.fetch_add(1, std::memory_order_acq_rel) + 1);
            //if((temp%10000) == 0) std::cout<<"Stored Steps: "<<temp<<std::endl;
            //check if we are done sampling (since we can only finish sampling on a written point
            if(stepsToTake == temp)
            {
                endOfSampling = true;
            }
            //check if we need to shift to skip mode for the next X points
            if(stepsToSkip != 0)
            {//whenever we write a point we need to skip the next batch before writing again
                savePoints.store(false, std::memory_order_release);
                stepsSkipped.store(0, std::memory_order_release);
            }
        }
        else
        {
            int temp = (stepsSkipped.fetch_add(1, std::memory_order_acq_rel) + 1);
            //we cannot be done sampling if we are in a stretch of skipping
            //check if we need to save the next step
            if(stepsToSkip == temp)
            {
                savePoints.store(false, std::memory_order_release);
            }
        }
        //handle if we are at the end of sampling
        if(endOfSampling)
        {
            numWorkersAtMainWait.store(0, std::memory_order_release);
            wrkrStatus.store(WorkerStatus::Wait, std::memory_order_release);
            std::unique_lock<std::mutex> ctrlLock(ctrlMutex);
            ctrlStatus.store(MainStatus::Continue, std::memory_order_release);
            ctrlWait.notify_all();
        }
        else
        {
            numWorkersAtMidStep = 0;
            wrkrIndex.store(0, std::memory_order_release);
            wrkrStatus.store(WorkerStatus::ProcessRed, std::memory_order_release);
        }
    }
    else
    {//we are not the last thread, spin until the state is not ProcessBlack
        while(wrkrStatus.load(std::memory_order_acquire) == WorkerStatus::ProcessBlack) {}
    }
}

template<class ParamType, class EndOfStepAction>
void RedBlackCtrler<ParamType, EndOfStepAction>::workerAckTerminate()
{
    int temp = wrkrTerm.fetch_add(1, std::memory_order_acq_rel);
    if(temp == lastThread)
    {
        std::unique_lock<std::mutex> ctrlLock(ctrlMutex);
        ctrlWait.notify_all();
    }
}

}
}
#endif  //MCMCPP_THREADING_REDBLKCTRLER_H
