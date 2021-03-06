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
#ifndef MCMCPP_CHAIN_CHAINPSETITERATOR_H
#define MCMCPP_CHAIN_CHAINPSETITERATOR_H
// includes for C system headers
// includes for C++ system headers
// includes from other libraries
// includes from MCMCpp
#include"ChainBlock.h"

namespace MCMC
{

namespace Chain
{

/**
 * @class ChainPsetIterator
 * @ingroup Chain
 * @brief In iterator to access a single parameter set from a single walker at a time
 * @author James Till Matta
 * 
 * \tparam ParamType The floating point type to be used for the chain, float, double, long double, etc.
 * 
 * ChainPsetIterator allows traversal point by point. That
 * is to say that incrementing the iterator moves to the next parameter set,
 * which could be in the same step just the next walker, or it could be in the
 * first walker of the next step.
 * 
 * ChainPsetIterator is a BidirectionalIterator, you can increment and decrement
 * it, as well as dereference it. However, you must use prefix increment and
 * decrement post-fix increment and decrement are not supported.
 * 
 * @remark It is possible to decrement an "end" iterator to yield an iterator
 * pointing to the very last parameter set in the last step taken
 */
template <class ParamType>
class ChainPsetIterator
{
public:
    /*!
     * \brief ChainPsetIterator Constructor
     * \param block The block the iterator is starting pointed to
     * \param blockStep The step within the block that the iterator starts pointing to
     */
    ChainPsetIterator(ChainBlock<ParamType>* block, int blockStep):
        curr(block), index(blockStep), lastCell(block->firstEmptyStep*block->walkerCount - 1),
        endCell(Detail::BlockSize*block->walkerCount - 1)
    {}
    
    /*!
     * \brief ChainPsetIterator Copy constructor to make a copy of an iterator
     * \param copy The original iterator to be copied into the iterator being constructed
     */
    ChainPsetIterator(const ChainPsetIterator<ParamType>& copy):
        curr(copy.curr), index(copy.index), lastCell(copy.lastCell),
        endCell(copy.endCell) {}
    
    ~ChainPsetIterator(){}
    
    /*!
     * \brief operator = Assignment operator for making two iterators that point to the same place
     * \param rhs The iterator to copy into this iterator
     * \return A reference to this iterator (post copying of rhs into it)
     */
    ChainPsetIterator<ParamType>& operator=(const ChainPsetIterator<ParamType>& rhs)
    {curr = rhs.curr; index = rhs.index; lastCell = rhs.lastCell; endCell = rhs.endCell; return *this;}
    
    /*!
     * \brief operator== equality operator to test equality of two iterators (equivalence of their locations)
     * \param rhs the second iterator (the right hand side)
     * \return  True if they are equal, false otherwise
     */
    bool operator==(const ChainPsetIterator<ParamType>& rhs){return ((curr == rhs.curr) && (index == rhs.index));}
    
    /*!
     * \brief operator!= equality operator to test inequality of two iterators (inequivalence of their locations)
     * \param rhs the second iterator (the right hand side)
     * \return  True if they are not equal, false otherwise
     */
    bool operator!=(const ChainPsetIterator<ParamType>& rhs){return ((curr != rhs.curr) || (index != rhs.index));}
    
    /*!
     * \brief operator++ Prefix increment of the iterator, move it to the next walker parameter set in the chain
     * \return The iterator that was incremented
     */
    ChainPsetIterator<ParamType>& operator++();
    /*!
     * \brief operator++ Prefix decrement of the iterator, move it to the previous walker parameter set in the chain
     * \return The iterator that was decremented
     */
    ChainPsetIterator<ParamType>& operator--();
    
    /*!
     * \brief operator* Dereference the iterator to get a pointer to the walker parameter array for this walker parameter set
     * \return A pointer into the walker parameter array, pointed to the beginning of the current walker's parameter set
     * 
     * @remark This function will let you dereference an invalid iterator, this
     * is undefined behaviour.
     */
    ParamType* operator*(){return (curr->chainArray + index*(curr->cellsPerWalker));}
    
    template<class> friend class Chain;
    
private:
    //Linked list book-keeping
    ChainBlock<ParamType>* curr = nullptr; ///<pointer to the current block
    int index;///<parameter set index within the current block
    int lastCell;///<stored the inddex of the last cell with data within the block array
    int endCell;///<stores the index of the last valid cell within the block array
};

template <class ParamType>
ChainPsetIterator<ParamType>& ChainPsetIterator<ParamType>::operator++()
{
    //check if we are not at the end of a block
    if(index < lastCell)
    {//not at the end of a block, easy peasy
        ++index;
        return *this;
    }
    else if((index>endCell) && (curr->nextBlock != nullptr))
    {//at the end of a block, there is a next block, jump to the next block
        curr = curr->nextBlock;
        index = 0;
        lastCell = (curr->firstEmptyStep*curr->walkerCount - 1);
        return *this;
    }
    else
    {//at the end of the chain
        index = (lastCell + 1);
        return *this;
    }
}

template <class ParamType>
ChainPsetIterator<ParamType>& ChainPsetIterator<ParamType>::operator--()
{
    //check if we are not at the end of a block
    if(index > 0)
    {//not at the start of a block, easy peasy
        --index;
        return *this;
    }
    else if(curr->nextBlock != nullptr)
    {//at the start of a block, jump to the end of the previous block
        curr = curr->lastBlock;
        lastCell = (curr->firstEmptyStep*curr->walkerCount - 1);
        index = lastCell;
        return *this;
    }
    else
    {//at the start of the chain
        index = 0;
        return *this;
    }
}

}
}
#endif  //MCMCPP_CHAIN_CHAINPSETITERATOR_H
