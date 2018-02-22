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
#ifndef MCMCPP_THREADING_REDBLKUPDATER_H
#define MCMCPP_THREADING_REDBLKUPDATER_H
// includes for C system headers
// includes for C++ system headers
#include<tuple>
// includes from other libraries
// includes from MCMCpp
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
template<class ParamType, class MoverType, class EndOfStepAction, bool UseSpinLocks>
class RedBlkUpdater
{
public:
    typedef std::tuple<Walker::Walker<ParamType>*, Walker::Walker<ParamType>*, int, int> WalkerInfo;
    typedef Walker::Walker<ParamType> WalkerType;
    typedef RedBlackCtrler<ParamType, EndOfStepAction, UseSpinLocks> ControllerType;
    
    /*!
     * \brief RedBlkUpdater Constructor for the Red Black updater thread
     * \param rndSeed The seed for the random number generator
     * \param threadNum The index of this thread, also the stream number for the rng
     * \param walkerSets A tuple of Walker::Walker<ParamType>*, Walker::Walker<ParamType>*, int, and int. The first pointer
     * is a pointer to the full set of red walkers, the second pointer is a pointer to the full set of black walkers, the first
     * integer is the total number of red walkers, the second integer is the total number of black walkers
     * \param mvr The mover class that will be copied to become the mover for this particular thread
     * \param ctrl A reference to the RedBlack controller (which will be stored)
     */
    RedBlkUpdater(int rndSeed, int threadNum, WalkerInfo& walkerSets, const MoverType& mvr, ControllerType& ctrl);
    
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
    int redSize = 0; ///<The number of walkers in the full red set
    int blkSize = 0; ///<The number of walkers in the full black set
    int trdNum; ///<The index of this thread, for debugging.
    bool writeChain = true;///<A variable that stores if the points need to be written for a given step
};

template<class ParamType, class MoverType, class EndOfStepAction, bool UseSpinLocks>
RedBlkUpdater<ParamType, MoverType, EndOfStepAction, UseSpinLocks>::
RedBlkUpdater(int rndSeed, int threadNum, WalkerInfo& walkerSets, const MoverType& mvr, ControllerType& ctrl):
    mover(mvr), controller(ctrl), trdNum(threadNum)
{
    std::tie(redWalkers, blkWalkers, redSize, blkSize) = walkerSets;
    mover.setPrng(rndSeed, threadNum);
}

template<class ParamType, class MoverType, class EndOfStepAction, bool UseSpinLocks>
void RedBlkUpdater<ParamType, MoverType, EndOfStepAction, UseSpinLocks>::operator()()
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
            notTerminated = false;
            break;
        }
    }
    //tell the controller that we have terminated
    controller.workerAckTerminate();
}

template<class ParamType, class MoverType, class EndOfStepAction, bool UseSpinLocks>
void RedBlkUpdater<ParamType, MoverType, EndOfStepAction, UseSpinLocks>::doRedUpdates()
{
    int index = controller.workerGetWorkIndex();
    //std::cout<<"RT"<<trdNum<<", "<<index<<std::endl;
    while(index<redSize)
    {
        mover.updateWalker(redWalkers[index], blkWalkers, blkSize, writeChain);
        index = controller.workerGetWorkIndex();
    }
}

template<class ParamType, class MoverType, class EndOfStepAction, bool UseSpinLocks>
void RedBlkUpdater<ParamType, MoverType, EndOfStepAction, UseSpinLocks>::doBlackUpdates()
{
    int index = controller.workerGetWorkIndex();
    //std::cout<<"BT"<<trdNum<<", "<<index<<std::endl;
    while(index < redSize)
    {
        mover.updateWalker(blkWalkers[index], redWalkers, redSize, writeChain);
        index = controller.workerGetWorkIndex();
    }
}

}
}
#endif  //MCMCPP_THREADING_REDBLKUPDATER_H
