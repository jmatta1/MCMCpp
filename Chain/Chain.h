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
#ifndef MCMC_CHAIN_CHAIN_H
#define MCMC_CHAIN_CHAIN_H
// includes for C system headers
// includes for C++ system headers
// includes from other libraries
// includes from MCMC
#include"ChainBlock.h"
#include"ChainPsetIterator.h"
#include"ChainStepIterator.h"

namespace MarkovChainMonteCarlo
{

namespace Chain
{

/*!
 * \brief The IncrementStatus enum represents the possible statuses after the
 * application of an end of step increment
 */
enum class IncrementStatus : char {NormalIncrement, ///< Standard increment, no moving to a new block or anything else
                                   NewBlock, ///< Either jumped to or allocated a new block for the chain
                                   EndOfChain ///< Reached the end of the last block *and* cannot allocate more memory to make a new block
                                  };

/**
 * @class Chain
 * @ingroup Chain
 * @brief The main access point of the chain which manages the linked list that underlies the chain.
 * @author James Till Matta
 * 
 * \tparam ParamType The floating point type to be used for the chain, float, double, long double, etc.
 * \tparam BlockSize The number of steps per walker that the block will hold
 * 
 * This class manages the chain of parameters. It allocates new nodes of the
 * linked list and passes insert parameters commands to the approporiate nodes.
 * It also can construct iterators that allow two styles of traversal of the
 * chain for later analysis. Usage of this class should be that of 'long term'
 * storage. That is, walkers maintain their own internal parameter sets, and
 * 'push' those sets to the chain as needed, as opposed to walkers attempting
 * to use memory allocated in the chain to store their current position
 * 
 * the first iterator ChainStepIterator allows traversal step by step. That is
 * to say that incrementing the iterator jumps all the parameter sets of all the
 * walkers, which is useful for operations that need to look at all the walkers
 * in a given point, at once.
 * 
 * The second iterator ChainPsetIterator allows traversal point by point. That
 * is to say that incrementing the iterator moves to the next parameter set,
 * which could be in the same step just the next walker, or it could be in the
 * first walker of the next step.
 */
template <class ParamType, int BlockSize>
class Chain
{
public:
    typedef ChainPsetIterator<ParamType, BlockSize> PsetIterator;
    typedef ChainStepIterator<ParamType, BlockSize> StepIterator;
    
    /*!
     * \brief Chain Constructor
     * \param numWalkers Number of walkers to store parameters for
     * \param numCellsPerWalker Number of parameters plus overhead that the walker needs
     * \param maxSize Maximum size in bytes that the chain is allowed to occupy (does not include overhead of nodes etc)
     */
    Chain(int numWalkers, int numCellsPerWalker, unsigned long long maxSize);
    ~Chain();
    
    /*!
     * \brief storeWalker Stores the walker's parameter set in the current step
     * \param walkerNum The index of the walker
     * \param walkerData The data array of the walker
     */
    void storeWalker(int walkerNum, ParamType* walkerData){curr->storeWalker(walkerNum, walkerData);}
    
    /*!
     * \brief incrementChainStep Moves the chain to the next step
     * \return An IncrementStatus enum whose value describes the chain state after the increment
     */
    IncrementStatus incrementChainStep();
    
    /*!
     * \brief getStoredStepCount get the number of steps
     * \return The number of steps stored in the chain
     */
    int getStoredStepCount(){return stepCount;}
    
    /*!
     * \brief incrementLastStep Moves the chain to the next step but will not allocate
     * a new block if the current step was the end of a block
     */
    void incrementLastStep(){++stepCount; curr->AndCheckNotFull();}
    
    /*!
     * \brief resetChain Effectively zeros out the number of steps taken by the walkers while leaving the memory allocated
     */
    void resetChain();
    
    /*!
     * \brief resetChainForSparse Resets the chain for sparse sampling post figuring out the auto correlation time
     * \param burnInSamples Number of samples at the start of the chain to discard in the name of 'burn-in'
     * \param autoCorrelationTime Longest autocorrelation time of all the parameters
     * 
     * This function discards the first numDiscardAc*autoCorrelationTime samples in the chain.
     * Of the remaining samples it then copies a sample every autoCorrelationTime samples along the chain
     * this places a set of completely independent samples of the distribution at the beginning of the 
     * chain and then leaves the remainder of the chain empty.
     */
    void resetChainForSubSampling(int burnInSamples, int autoCorrelationTime);
    
    /*!
     * \brief getPsetIteratorBegin Gets a parameter set iterator pointed at the very first parameter set
     * \return ChainPsetIterator pointed to the very beginning of the chain's parameter sets
     */
    ChainPsetIterator<ParamType, BlockSize> getPsetIteratorBegin(){return PsetIterator(head, 0);}
    /*!
     * \brief getPsetIteratorEnd Gets a parameter set iterator pointed just after the very last parameter set
     * \return ChainPsetIterator pointed to just past the last of the chain's parameter sets
     */
    ChainPsetIterator<ParamType, BlockSize> getPsetIteratorEnd();
    
    /*!
     * \brief getStepIteratorBegin Gets a step iterator pointed at the very first step in the chain
     * \return ChainStepIterator pointed to the beginning of the chain
     */
    ChainStepIterator<ParamType, BlockSize> getStepIteratorBegin(){return StepIterator(head, 0);}
    /*!
     * \brief getStepIteratorEnd Gets a step iterator pointed to just past the end of the chain
     * \return ChainStepIterator pointed to just after the end of the chain
     */
    ChainStepIterator<ParamType, BlockSize> getStepIteratorEnd();
private:
    /*!
     * \brief incrementChainStepSubSampleReset Moves the chain to the next step, doing special actions needed for resetting the chain as we go along
     */
    void incrementChainStepSubSampleReset();
    
    //Linked list book-keeping
    ChainBlock<ParamType, BlockSize>* head = nullptr; ///<pointer to the first block in the chain linked list
    /*!
     * \brief pointer to the last block in the chain linked list
     * 
     * if the chain is active, this will always point to a block that has a free
     * step in it, if the chain was finalized with incrementLastStep(), this
     * holds the final block in the chain, regardless of if there is a free
     * step in it or not
     */
    ChainBlock<ParamType, BlockSize>* curr = nullptr;
    //Chain book-keeping
    int walkerCount; ///<Number of walkers included in this chain
    int cellsPerWalker; ///<Number of cells needed by each walker
    int maxBlocks; ///<Calculated maximum number of blocks to allocate
    int stepCount = 0; ///<Number of steps stored so far
    int blockCount = 0; ///<Number of blocks allocated so far
};


template <class ParamType, int BlockSize>
Chain<ParamType, BlockSize>::Chain(int numWalkers, int numCellsPerWalker, unsigned long long maxSize):
    walkerCount(numWalkers), cellsPerWalker(numCellsPerWalker),
    maxBlocks(maxSize/(sizeof(ParamType)*BlockSize*numWalkers*cellsPerWalker))
{
    //allocate the first block
    head = new ChainBlock<ParamType, BlockSize>(nullptr, walkerCount, cellsPerWalker);
    //record that we have allocated a new block
    ++blockCount;
    //set up the curr ptr
    curr = head;
}

template <class ParamType, int BlockSize>
Chain<ParamType, BlockSize>::~Chain()
{
    //grab a pointer to head's next block, head should *always be valid
    ChainBlock<ParamType, BlockSize>* temp = head->nextBlock;
    delete head;
    // while temp is not a nullptr, move head to temp, move temp to nextblock and delete head, rinse and repeat
    while(temp != nullptr)
    {
        head = temp;
        temp = head->nextBlock;
        delete head;
    }
}

template <class ParamType, int BlockSize>
IncrementStatus Chain<ParamType, BlockSize>::incrementChainStep()
{
    ++stepCount;
    if(curr->incrementChainStepAndCheckNotFull())
    {//The increment when completely normally
        return IncrementStatus::NormalIncrement;
    }
    else
    {//The current block is now full we either need to go to the next block, allocate a new block, or declare the chain finished
        if(curr->nextBlock != nullptr)
        {//There is a next block already allocated, move to it
            curr = curr->nextBlock;
            //Tell the caller that the chain has moved to a new block
            return IncrementStatus::NewBlock;
        }
        else if(blockCount < maxBlocks)
        {//There is no pre-allocated next block, but we can allocate one
            curr->nextBlock = new ChainBlock<ParamType, BlockSize>(curr, walkerCount, cellsPerWalker);
            //record the additional allocation
            ++blockCount;
            //move to the new block
            curr = curr->nextBlock;
            //tell the caller that the chain has moved to a new block
            return IncrementStatus::NewBlock;
        }
        else
        {//there is no next block available and we cannot allocate anymore without exceeding the space limits
            //tell the user
            return IncrementStatus::EndOfChain;
        }
    }
}

template <class ParamType, int BlockSize>
void Chain<ParamType, BlockSize>::incrementChainStepSubSampleReset()
{
    //increment the number of steps taken
    ++stepCount;
    if(!curr->incrementChainStepAndCheckNotFull())
    {
        //The current block is now full we either need to go to the next block
        if(curr->nextBlock != nullptr)
        {//There is a next block already allocated, move to it
            curr = curr->nextBlock;
            //since this function is only called by resetChainForSubSampling we need to reset the new block we moved into
            curr->reset();
        }//The else cases should *NOT* be able to happen when called from resetChainForSubSampling
    }
}

template <class ParamType, int BlockSize>
void Chain<ParamType, BlockSize>::resetChain()
{
    ChainBlock<ParamType, BlockSize>* temp = head;
    while(temp != nullptr)
    {
        temp->reset();
        temp = temp->nextBlock;
    }
    stepCount = 0;
}

template <class ParamType, int BlockSize>
void Chain<ParamType, BlockSize>::resetChainForSubSampling(int burnInSamples, int autoCorrelationTime)
{
    stepCount = 0;
    auto readLocation = this->getStepIteratorBegin();
    auto end = this->getStepIteratorBegin();
    //push the read location to the first non-burnin sample
    readLocation += burnInSamples;
    //check to make sure we are not at the end (if so, just do a normal reset of the chain
    if(readLocation == end){resetChain(); return;}
    //if we are here we are not at the end of the chain and so we need to start copying samples properly
    //first put our current pointer at head
    curr = head;
    //now reset head
    //it is safe to reset head, even though the iterator needs to know the size (incase the iterator is within head still)
    //this is because the iterator gets what it needs as it is constructed. Post Construction it updates itself when it jumps to
    //the new block, since the point being written will always be behind the point being read, the iterator will always get what
    //it needs before the block is reset
    curr->reset();
    //now loop and copy walker parameter sets
    while(reeadLocation != end)
    {
        curr->copyWalkerSet(*readLocation);
        incrementChainStepSubSampleReset();
        readLocation += autoCorrelationTime;
    }
    //now finish resetting the remainder of the chain
    ChainBlock<ParamType, BlockSize>* temp = curr->nextBlock;
    while(temp != nullptr)
    {
        temp->reset();
        temp = temp->nextBlock;
    }
}

template <class ParamType, int BlockSize>
ChainPsetIterator<ParamType, BlockSize> Chain<ParamType, BlockSize>::getPsetIteratorEnd()
{
    int index = (curr->firstEmptyStep*walkerCount);
    return ChainPsetIterator<ParamType, BlockSize>(curr, index);
}

template <class ParamType, int BlockSize>
ChainStepIterator<ParamType, BlockSize> Chain<ParamType, BlockSize>::getStepIteratorEnd()
{
    return ChainPsetIterator<ParamType, BlockSize>(curr, curr->firstEmptyStep);
}

}
}
#endif  //MCMC_CHAIN_CHAIN_H
