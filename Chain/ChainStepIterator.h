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
 * ChainStepIterator is a BidirectionalIterator, you can increment and decrement
 * it, as well as dereference it. However, you must use prefix increment and
 * decrement, post-fix increment and decrement are not supported.
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
    ~ChainStepIterator(){}
    
    /*!
     * \brief operator== equality operator to test equality of two iterators (equivalence of their locations)
     * \param rhs the second iterator (the right hand side)
     * \return  True if they are equal, false otherwise
     */
    bool operator==(const ChainStepIterator& rhs){return ((curr == rhs.curr) && (stepIndex == rhs.stepIndex));}
    /*!
     * \brief operator!= equality operator to test inequality of two iterators (inequivalence of their locations)
     * \param rhs the second iterator (the right hand side)
     * \return  True if they are not equal, false otherwise
     */
    bool operator!=(const ChainStepIterator& rhs){return ((curr != rhs.curr) || (stepIndex != rhs.stepIndex));}
    
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
    int stepIndex;///<step index within the current block
};

template <class ParamType, int BlockSize>
ChainStepIterator<ParamType, BlockSize> ChainStepIterator<ParamType, BlockSize>::operator++()
{
    /*!
     * @remark This function will allow the iterator to be incremented past the
     * end. In fact, if the end is not the last cell of the last block merely a
     * a cell in the last block, then it is possible to increment many times past
     * the end iterator, once the iterator has gone past the last cell of the
     * last block, the iterator is not recoverable
     */
    //check if we are not at the end of a block
    if(stepIndex < (BlockSize-1))
    {//not at the end of a block, easy peasy
        ++stepIndex;
    }
    else
    {//at the end of a block, jump to the next block
        curr = curr->nextBlock;
        stepIndex = 0;
    }
    return *this;
}

template <class ParamType, int BlockSize>
ChainStepIterator<ParamType, BlockSize> ChainStepIterator<ParamType, BlockSize>::operator--()
{
    /*!
     * @remark This function will allow a .begin() operator to be decremented,
     * once this happens the iterator is not recoverable
     */
    //check if we are not at the end of a block
    if(stepIndex > 0)
    {//not at the start of a block, easy peasy
        --stepIndex;
    }
    else
    {//at the start of a block, jump to the end of the previous block
        curr = curr->lastBlock;
        stepIndex = (BlockSize-1);
    }
    return *this;
}

}
}
#endif  //MCMC_CHAIN_CHAINSTEPITERATOR_H
