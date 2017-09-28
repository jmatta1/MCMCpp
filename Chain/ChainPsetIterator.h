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
#ifndef MCMC_CHAIN_CHAINPSETITERATOR_H
#define MCMC_CHAIN_CHAINPSETITERATOR_H
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
 * @brief In iterator to access all the parameter sets of a group of walkers simultaneously
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
 * ChainStep is a BidirectionalIterator, you can increment and decrement it, as
 * well as dereference it. However, you must use prefix increment and decrement
 * post-fix increment and decrement are not supported.
 */
template <class ParamType, int BlockSize>
class ChainPsetIterator
{
public:
    /*!
     * \brief Chain Constructor
     * \param block The block the iterator is starting pointed to
     * \param blockStep The step within the block that the iterator starts pointing to
     */
    ChainPsetIterator(ChainBlock<ParamType, BlockSize>* block, int blockStep):
        curr(block), stepIndex(blockStep), lastCell(BlockSize*block->walkerCount - 1){}
    ~ChainPsetIterator(){}
    
    bool operator==(const ChainStepIterator& rhs){return ((curr == rhs.curr) && (index == rhs.index));}
    bool operator!=(const ChainStepIterator& rhs){return ((curr != rhs.curr) || (index != rhs.index));}
    
    ChainPsetIterator<ParamType, BlockSize> operator++();
    ChainPsetIterator<ParamType, BlockSize> operator--();
    
    ParamType* operator*(){return (curr->chainArray + index*(curr->cellsPerWalker));}
    
private:
    //Linked list book-keeping
    ChainBlock<ParamType, BlockSize>* curr = nullptr; ///<pointer to the current block
    int index;///<parameter set index within the current block
    int lastCell;///<stores the index of the last valid cell within the block array
};

template <class ParamType, int BlockSize>
ChainPsetIterator<ParamType, BlockSize> ChainPsetIterator<ParamType, BlockSize>::operator++()
{
    //check if we are not at the end of a block
    if(index < lastCell)
    {//not at the end of a block, easy peasy
        ++index;
    }
    else
    {//at the end of a block, jump to the next block
        curr = curr->nextBlock;
        index = 0;
    }
    return *this;
}

template <class ParamType, int BlockSize>
ChainPsetIterator<ParamType, BlockSize> ChainPsetIterator<ParamType, BlockSize>::operator--()
{
    //check if we are not at the end of a block
    if(index > 0)
    {//not at the start of a block, easy peasy
        --index;
    }
    else
    {//at the start of a block, jump to the end of the previous block
        curr = curr->lastBlock;
        index = lastCell;
    }
    return *this;
}

}
}
#endif  //MCMC_CHAIN_CHAINPSETITERATOR_H
