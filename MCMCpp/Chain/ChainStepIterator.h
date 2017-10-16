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
#ifndef MCMC_CHAIN_CHAINSTEPITERATOR_H
#define MCMC_CHAIN_CHAINSTEPITERATOR_H
// includes for C system headers
// includes for C++ system headers
// includes from other libraries
// includes from MCMC
#include"ChainBlock.h"

namespace MarkovChainMonteCarlo
{

namespace Chain
{

/**
 * @class ChainStepIterator
 * @ingroup Chain
 * @brief In iterator to access all the parameter sets of the set of walkers simultaneously
 * @author James Till Matta
 * 
 * \tparam ParamType The floating point type to be used for the chain, float, double, long double, etc.
 * \tparam BlockSize The number of steps per walker that the block will hold
 * 
 * ChainStepIterator allows traversal of the chain step by step. That is
 * to say that incrementing the iterator jumps all the parameter sets of all the
 * walkers, which is useful for operations that need to look at all the walkers
 * in a given point, at once.
 * 
 * ChainStepIterator is a bit more than a BidirectionalIterator, you can
 * increment (prefix only) and decrement (prefix only) it, as well as
 * dereference it. However, you can also add and subtract arbitrary numbers of
 * steps with += and -=. It is not a RandomAccess iterator for a few reasons
 * though. The += and -= are not fixed time and  there is no support for the
 * following operators >, >=, <, <=, -, +, or []. The += and -= are provided
 * because they allow much more efficient jumps than simply using ++ or --
 * n times.
 * 
 * @remark It is possible to decrement an "end" iterator to yield an iterator
 * pointing to the very last step taken
 */
template <class ParamType, int BlockSize>
class ChainStepIterator
{
public:
    /*!
     * \brief ChainStepIterator Constructor
     * \param block The block the iterator is starting pointed to
     * \param blockStep The step within the block that the iterator starts pointing to
     */
    ChainStepIterator(ChainBlock<ParamType, BlockSize>* block, int blockStep):curr(block), stepIndex(blockStep){}
    /*!
     * \brief ChainStepIterator Copy constructor to make a copy of an iterator
     * \param copy The original iterator to be copied into the iterator being constructed
     */
    ChainStepIterator(const ChainStepIterator<ParamType, BlockSize>& copy):
        curr(copy.curr), lastFullStep(copy.lastFullStep), stepIndex(copy.stepIndex)
    {}
    
    ~ChainStepIterator(){}
    /*!
     * \brief operator = Assignment operator for making two iterators that point to the same place
     * \param rhs The iterator to copy into this iterator
     * \return A reference to this iterator (post copying of rhs into it)
     */
    ChainStepIterator<ParamType, BlockSize>& operator=(const ChainStepIterator<ParamType, BlockSize>& rhs)
    {curr = rhs.curr; lastFullStep = rhs.lastFullStep; stepIndex = rhs.stepIndex; return *this;}
    
    /*!
     * \brief operator== equality operator to test equality of two iterators (equivalence of their locations)
     * \param rhs the second iterator (the right hand side)
     * \return  True if they are equal, false otherwise
     */
    bool operator==(const ChainStepIterator<ParamType, BlockSize>& rhs){return ((curr == rhs.curr) && (stepIndex == rhs.stepIndex));}
    /*!
     * \brief operator!= equality operator to test inequality of two iterators (inequivalence of their locations)
     * \param rhs the second iterator (the right hand side)
     * \return  True if they are not equal, false otherwise
     */
    bool operator!=(const ChainStepIterator<ParamType, BlockSize>& rhs){return ((curr != rhs.curr) || (stepIndex != rhs.stepIndex));}
    
    /*!
     * \brief operator++ Prefix increment of the iterator, move it to the next step in the chain
     * \return The iterator that was incremented
     */
    ChainStepIterator<ParamType, BlockSize> operator++();
    /*!
     * \brief operator++ Prefix decrement of the iterator, move it to the previous step in the chain
     * \return The iterator that was decremented
     */
    ChainStepIterator<ParamType, BlockSize> operator--();
    
    /*!
     * \brief operator+= Increase the iterator by some number of steps (stopping at the end if needed)
     * \return The iterator that was increased
     */
    ChainStepIterator<ParamType, BlockSize> operator+=(int steps);
    
    /*!
     * \brief operator-= Decrease the iterator by some number of steps (stopping at the Beginning if needed)
     * \return The iterator that was decreased
     */
    ChainStepIterator<ParamType, BlockSize> operator+=(int steps);
    
    /*!
     * \brief operator* Dereference the iterator to get a pointer to the walker parameter array for this step
     * \return A pointer into the walker parameter array, pointed to the beginning of this step
     * 
     * @remark This function will let you dereference an invalid iterator, this
     * is undefined behaviour.
     */
    ParamType* operator*(){return (curr->chainArray + stepIndex*(curr->cellsPerStep));}
    
    friend class Chain;
    
private:
    //Linked list book-keeping
    ChainBlock<ParamType, BlockSize>* curr = nullptr; ///<pointer to the current block
    int lastFullStep; ///<Proportional to the firstEmptyStep parameter of the current block
    int stepIndex;///<step index within the current block
};

template <class ParamType, int BlockSize>
ChainStepIterator<ParamType, BlockSize> ChainStepIterator<ParamType, BlockSize>::operator+=(int steps)
{
    /*!
     * @remark This increment *will* stop at the end of the chain, even if steps
     * would take it further than that, this is necessary for knowing when we
     * are at the end of the iterator
     */
    //first figure out how many blocks we need to jump
    int blocks = (steps/BlockSize);
    while((blocks > 0) && (curr->nextBlock != nullptr))
    {
        curr = curr->nextBlock;
        //leave step index alone because we are essentially incrementing by 1000
        lastFullStep = (curr->firstEmptyStep - 1);
        --blocks;
    }
    //Figure out if we stopped because we were at the end or because we hit our block Count
    if(blocks > 0)
    {//we stopped because we hit the end
        stepIndex = (lastFullStep + 1);
        return *this;
    }
    //if we are here we stopped because we hit our block count, now see if we need to jump one more block
    steps = (steps%BlockSize);
    if((stepIndex + steps) > BlockSize)
    {//we need to jump one more block
        if(curr->nextBlock != nullptr)
        {//there is a next block available
            curr = curr->nextBlock;
            //now we are incrementing by the remaining size of the last block,
            //therefor set step index to 0
            stepIndex = 0;
            lastFullStep = (curr->firstEmptyStep - 1);
            steps -= (BlockSize - stepIndex);
        }
        else
        {//the next block is not available
            stepIndex = (lastFullStep + 1);
            return *this;
        }
    }
    //if we are here, we do not need to jump another block, just add some steps
    //to the current step index
    stepIndex += steps;
    if(stepIndex <= lastFullStep)
    {
        return *this;
    }
    else
    {//we could possibly be in the final block of the chain and the steps is only a little more than is in the chain
        //if we are here then that is the case
        stepIndex = lastFullStep;
        return *this;
    }
}

template <class ParamType, int BlockSize>
ChainStepIterator<ParamType, BlockSize> ChainStepIterator<ParamType, BlockSize>::operator-=(int steps)
{
    /*!
     * @remark This increment *will* stop at the beginning of the chain, even if steps
     * would take it further than that, this is necessary for knowing when we
     * are at the end of the iterator
     */
    //first figure out how many blocks we need to jump
    int blocks = (steps/BlockSize);
    while((blocks > 0) && (curr->lastBlock != nullptr))
    {
        curr = curr->lastBlock;
        //leave step index alone because we are essentially incrementing by BlockSize
        lastFullStep = (curr->firstEmptyStep - 1);
        --blocks;
    }
    //Figure out if we stopped because we were at the end or because we hit our block Count
    if(blocks > 0)
    {//we stopped because we hit the end
        stepIndex = 0;
        return *this;
    }
    //if we are here we stopped because we hit our block count, now see if we need to jump one more block
    steps = (steps%BlockSize);
    if(stepIndex < steps)
    {//we need to jump one more block
        if(curr->lastBlock != nullptr)
        {//there is a last block available
            curr = curr->lastBlock;
            //now we are decrementing by the remaining size of the previous block,
            //therefor set step index to the last step in the new block
            stepIndex = (BlockSize-1);
            lastFullStep = (curr->firstEmptyStep - 1);
            steps -= stepIndex;
        }
        else
        {//We are in the first block
            stepIndex = 0;
            return *this;
        }
    }
    //if we are here, we do not need to jump another block, just add some steps
    //to the current step index
    stepIndex -= steps;
    if(stepIndex >= 0)
    {
        return *this;
    }
    else
    {//If we were in the first block of the chain and had enough steps to go past
        //the beginning we should have gotten stopped earlier in this function
        //this is here for completeness/symmetry with +=
        stepIndex = 0;
        return *this;
    }
}


template <class ParamType, int BlockSize>
ChainStepIterator<ParamType, BlockSize> ChainStepIterator<ParamType, BlockSize>::operator++()
{
    //check if we are not at the end of a block
    if(stepIndex < lastFullStep)
    {//not at the end of a block, easy peasy
        ++stepIndex;
    }
    else if((stepIndex==BlockSize) && (curr->nextBlock != nullptr))
    {//at the end of a block, and there is a next block, jump to the next block
        curr = curr->nextBlock;
        stepIndex = 0;
        lastFullStep = (curr->firstEmptyStep - 1);
    }
    else
    {//at the end of a block, and there is no next block, set step index to the first available step in the block
        stepIndex = (lastFullStep + 1);
    }
    return *this;
}

template <class ParamType, int BlockSize>
ChainStepIterator<ParamType, BlockSize> ChainStepIterator<ParamType, BlockSize>::operator--()
{
    //check if we are not at the end of a block
    if(stepIndex > 0)
    {//not at the start of a block, easy peasy
        --stepIndex;
    }
    else if(curr->lastBlock != nullptr)
    {//at the start of a block, and there is a previous block jump to the end of the previous block
        curr = curr->lastBlock;
        lastFullStep = (curr->firstEmptyStep - 1);
        stepIndex = lastFullStep;
    }
    else
    {//otherwise we are at the beginning of the first block
        stepIndex = 0;
    }
    return *this;
}

}
}
#endif  //MCMC_CHAIN_CHAINSTEPITERATOR_H
