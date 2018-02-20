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
#ifndef MCMCPP_CHAIN_CHAINBLOCK_H
#define MCMCPP_CHAIN_CHAINBLOCK_H
// includes for C system headers
// includes for C++ system headers
#include<cstdlib>//needed for aligned allocation
#include<cstring>//for memcpy
#include<algorithm>
// includes from other libraries
// includes from MCMCpp
#include"../Utility/Misc.h"

namespace MCMC
{

namespace Chain
{

namespace Detail
{
static const int BlockSize = 10000; ///<Number of steps to store in one chain block
}

/**
 * @class ChainBlock
 * @ingroup Chain
 * @brief Nodes of the linked list that represents the chain
 * @author James Till Matta
 * 
 * \tparam ParamType The floating point type to be used for the chain, float, double, long double, etc.
 * 
 * This class holds 'chunks' of the chain that is distributed across a linked
 * list style structure. Each block holds Detail::BlockSize parameter sets for all the
 * walkers
 * 
 */
template <class ParamType>
class ChainBlock
{
public:
    /*!
     * \brief ChainBlock Constructor
     * \param prev Previous block in the chain, can be nullptr for the head node
     * \param numWalkers Number of walkers to store parameters for
     * \param numCellsPerWalker Number of parameters plus overhead that the walker needs
     */
    ChainBlock(ChainBlock<ParamType>* prev, int numWalkers, int numCellsPerWalker);
    ~ChainBlock(){Utility::delAAA(chainArray);}
    
    /*!
     * \brief Deleted copy constructor
     */
    ChainBlock(const ChainBlock<ParamType>& rhs) = delete;
    
    /*!
     * \brief Deleted assignment operator
     */
    ChainBlock<ParamType>& operator=(const ChainBlock<ParamType>& rhs) = delete;
    
    /*!
     * \brief storeWalker transfers a walkers point into the chain storage, undefined behavior if the block is full
     * \param walkerNum index of the walker to transfer
     * \param walkerData array of data representing the walker's parameters and auxilliary data
     * 
     * @remark This function should be reentrant, that is, two simultaneous calles with different walkerNum will not have problems
     */
    void storeWalker(int walkerNum, ParamType* walkerData);
    
    /*!
     * \brief incrementChainStepAndCheckSpace should be called after a the full set of walkers has stored themselves to the chain
     * \return Returns true if there is space for at least one more step in the block, otherwise returns false
     */
    bool incrementChainStepAndCheckNotFull(){++firstEmptyStep; return (firstEmptyStep < Detail::BlockSize);}
    
    /*!
     * \brief copyWalkerSet Copys a full set of walker parameters from the given location into the current step
     * \param walkerData Array starting at the beginning of a step for a full set of walkers
     */
    void copyWalkerSet(ParamType* walkerData);
    
    /*!
     * \brief reset Resets the block so it contains no steps
     */
    void reset(){firstEmptyStep = 0;}
    
    template<class> friend class Chain;
    
    template<class> friend class ChainStepIterator;
    
    template<class> friend class ChainPsetIterator;
private:
    //Linked list book-keeping
    ChainBlock<ParamType>* lastBlock; ///<pointer to the previous block in the chain linked list
    ChainBlock<ParamType>* nextBlock = nullptr; ///<pointer to the next block in the chain linked list
    //the chain chunk
    ParamType* chainArray;///<Pointer into the array that stores the parameter set for every walker for every step in the block
    //Chain book-keeping
    int walkerCount; ///<Number of walkers included in this chain
    int cellsPerWalker; ///<Number of cells needed by each walker
    int cellsPerStep; ///<Number of cells for all the walkers in a single step
    int firstEmptyStep = 0; ///<Index (in units of cellsPerStep) of the next empty step
    int sizeInBytes; ///<The size in bytes for a memcpy operation
};

template <class ParamType>
ChainBlock<ParamType>::ChainBlock(ChainBlock<ParamType>* prev, int numWalkers, int numCellsPerWalker):
    lastBlock(prev),chainArray(nullptr),
    walkerCount(numWalkers), cellsPerWalker(numCellsPerWalker),
    cellsPerStep(numWalkers*numCellsPerWalker),
    sizeInBytes(walkerCount*cellsPerWalker*sizeof(ParamType))
{
    chainArray = Utility::autoAlignedAlloc<ParamType>(Detail::BlockSize*numCellsPerWalker*numWalkers);
}

template <class ParamType>
void ChainBlock<ParamType>::storeWalker(int walkerNum, ParamType* walkerData)
{
    ParamType* offset = chainArray + (firstEmptyStep*cellsPerStep + walkerNum*cellsPerWalker);
    //std::copy(walkerData, (walkerData+cellsPerWalker), offset);
    std::memcpy(offset, walkerData, cellsPerWalker*sizeof(ParamType));
}

template <class ParamType>
void ChainBlock<ParamType>::copyWalkerSet(ParamType* walkerData)
{
    ParamType* offset = chainArray + (firstEmptyStep*cellsPerStep);
    //std::copy(walkerData, (walkerData+(walkerCount*cellsPerWalker)), offset);
    std::memcpy(offset, walkerData, sizeInBytes);
}

}
}
#endif  //MCMCPP_CHAIN_CHAINBLOCK_H
