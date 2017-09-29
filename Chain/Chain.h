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
     * \return True if there is another step available (i.e. maxBlock has been reached), False otherwise
     */
    bool incrementChainStep();
    
    /*!
     * \brief getStoredStepCount get the number of steps
     * \return The number of steps stored in the chain
     */
    int getStoredStepCount(){return stepCount;}
    
    /*!
     * \brief incrementLastStep Moves the chain to the next step but will not allocate
     * a new block if the current step was the end of a block
     */
    void incrementLastStep(){++stepCount; curr->incrementChainStepAndCheckNotFull();}
    
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
    ++blockCount;
    curr = head;
}

template <class ParamType, int BlockSize>
Chain<ParamType, BlockSize>::~Chain()
{
    ChainBlock<ParamType, BlockSize>* temp = head->nextBlock;
    delete head;
    while(temp != nullptr)
    {
        head = temp;
        temp = head->nextBlock;
        delete head;
    }
}

template <class ParamType, int BlockSize>
bool Chain<ParamType, BlockSize>::incrementChainStep()
{
    ++stepCount;
    if(!curr->incrementChainStepAndCheckNotFull())
    {
        if(blockCount < maxBlocks)
        {
            curr->nextBlock = new ChainBlock<ParamType, BlockSize>(curr, walkerCount, cellsPerWalker);
            ++blockCount;
            curr = curr->nextBlock;
            return true;
        }
        else
        {
            return false;
        }
    }
    else
    {
        return true;
    }
}

template <class ParamType, int BlockSize>
ChainPsetIterator<ParamType, BlockSize> Chain<ParamType, BlockSize>::getPsetIteratorEnd()
{
    if(curr->firstEmptyStep == BlockSize)
    {
        return ChainPsetIterator<ParamType, BlockSize>(nullptr, 0);
    }
    else
    {
        int index = (curr->firstEmptyStep*walkerCount*cellsPerWalker);
        return ChainPsetIterator<ParamType, BlockSize>(curr, index);
    }
}

template <class ParamType, int BlockSize>
ChainStepIterator<ParamType, BlockSize> Chain<ParamType, BlockSize>::getStepIteratorEnd()
{
    if(curr->firstEmptyStep == BlockSize)
    {
        return ChainPsetIterator<ParamType, BlockSize>(nullptr, 0);
    }
    else
    {
        return ChainPsetIterator<ParamType, BlockSize>(curr, curr->firstEmptyStep);
    }
}

}
}
#endif  //MCMC_CHAIN_CHAIN_H
