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
#ifndef MCMC_THREADING_REDBLKUPDATER_H
#define MCMC_THREADING_REDBLKUPDATER_H
// includes for C system headers
// includes for C++ system headers
#include<tuple>
// includes from other libraries
// includes from MCMC
#include"../Walker/Walker.h"
#include"RedBlkCtrler.h"

namespace MCMC
{
namespace Threading
{

/**
 * @class RedBlackUpdater
 * @ingroup Threading
 * @brief The callable object used to run an updater thread
 * @author James Till Matta
 * 
 */
template<class ParamType, class MoverType, class EndOfStepAction>
class RedBlkUpdater
{
public:
    typedef std::tuple<Walker::Walker<ParamType>*, Walker::Walker<ParamType>*, int, int> WalkerInfo;
    typedef Walker::Walker<ParamType> WalkerType;
    typedef RedBlackCtrler<ParamType, EndOfStepAction> ControllerType;
    
    /*!
     * \brief RedBlkUpdater Constructor for the Red Black updater thread
     * \param rndSeed The seed for the random number generator
     * \param threadNum The index of this thread, also the stream number for the rng
     * \param walkerSets A tuple of Walker::Walker<ParamType>*, Walker::Walker<ParamType>*, int, and int. The first pointer
     * is a pointer to the full set of red walkers, the second pointer is a pointer to the full set of black walkers, the first
     * integer is the total number of red walkers, the second integer is the total number of black walkers
     * \param updateSetsA tuple of Walker::Walker<ParamType>*, Walker::Walker<ParamType>*, int, and int. The first pointer
     * is a pointer to the set of red walkers that this thread is responsible for updating, the second pointer is a pointer
     * to the set of black walkers that this thread is responsible for updating. The first integer is the number of red
     * walkers that this thread is responsible for updating, the second integer is the number of black walkers that this
     * thread is responsible for updating
     * \param mvr The mover class that will be copied to become the mover for this particular thread
     * \param ctrl A reference to the RedBlack controller (which will be stored)
     */
    RedBlkUpdater(int rndSeed, int threadNum, WalkerInfo& walkerSets, WalkerInfo& updateSets, MoverType& mvr, ControllerType& ctrl);
    
    /*!
     * \brief operator() The entry point for the thread, carries out all the threads basic tasks
     */
    void operator()();
    
private:
    /*!
     * \brief doRedUpdates Updates the set of red walkers this thread is responsible for using the black walkers
     */
    void doRedUpdates();
    
    /*!
     * \brief doBlackUpdates Updates the set of black walkers that this thread is responsible for using the red walkers
     */
    void doBlackUpdates();
    
    MoverType mover; ///<A copy of the mover that will be used by this thread to process events
    ControllerType& controller; ///<A reference to the controller object
    WalkerType* redWalkers = nullptr; ///<A pointer to the beginning of the full set of red walkers, which will be used to update the black walkers
    WalkerType* blkWalkers = nullptr; ///<A pointer to the beginning of the full set of black walkers, which will be used to update the red walkers
    WalkerType* redSet = nullptr; ///<A pointer to the set of walkers that are updated for the red step
    WalkerType* blkSet = nullptr; ///<A pointer to the set of walkers that are updated for the black step
    int redSize = 0; ///<The number of walkers in the full red set
    int blkSize = 0; ///<The number of walkers in the full black set
    int numRed = 0; ///<The number of walkers in the red step that this thread is in charge of updating
    int numBlk = 0; ///<The number of walkers in the black step that this thread is in charge of updating
    
    bool writeChain = true;
};

template<class ParamType, class MoverType, class EndOfStepAction>
RedBlkUpdater<ParamType, MoverType, EndOfStepAction>::RedBlkUpdater(int runNum, int threadNum, WalkerInfo& walkerSets, WalkerInfo& updateSets, MoverType& mvr, ControllerType& ctrl):
    mover(mvr), controller(ctrl)
{
    std::tie(redWalkers, blkWalkers, redSize, blkSize) = walkerSets;
    std::tie(redSet, blkSet, numRed, numBlk) = updateSets;
    mover.setPrng(rndSeed, threadNum);
}

template<class ParamType, class MoverType, class EndOfStepAction>
void RedBlkUpdater<ParamType, MoverType, EndOfStepAction>::operator()()
{
    bool notTerminated = true;
    while(notTerminated)
    {
        WorkerStatus mode = controller.workerGetStatus();
        switch(mode)
        {
        case WorkerStatus::ProcessRed: //Red processing always comes first
            //retrieve if we are writing this point or not, we do not do this
            //for the black processing state because this state always came first
            writeChain = controller.workerGetSavingPoint();
            //do the updates for the red set walkers this thread is in charge of
            doRedUpdates();
            //now get the controller to figure out waiting for the next step / state
            controller.workerMidStepWait();
            break;
        case WorkerStatus::ProcessBlack:
            //do the updates for the black set walkers this thread is in charge of
            doBlackUpdates();
            //now get the controller to figure out waiting for the next step / state
            controller.workerEndStepWait();
            break;
        case WorkerStatus::Wait:
            controller.workerWaitForTask();
            break;
        case WorkerStatus::Terminate:
            notTerminated = true;
            break;
        }
    }
    //tell the controller that we have terminated
    controller.workerAckTerminate();
}

template<class ParamType, class MoverType, class EndOfStepAction>
void RedBlkUpdater<ParamType, MoverType, EndOfStepAction>::doRedUpdates()
{
    for(int i=0; i<numRed; ++i)
    {
        mover.updateWalker(redSet[i], blkWalkers, blkSize, writeChain);
    }
}

template<class ParamType, class MoverType, class EndOfStepAction>
void RedBlkUpdater<ParamType, MoverType, EndOfStepAction>::doBlackUpdates()
{
    for(int i=0; i<numBlk; ++i)
    {
        mover.updateWalker(blkSet[i], redWalkers, redSize, writeChain);
    }
}

}
}
#endif  //MCMC_THREADING_REDBLKUPDATER_H
