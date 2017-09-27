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
#ifndef MCMC_CHAIN_CHAINBLOCK_H
#define MCMC_CHAIN_CHAINBLOCK_H
// includes for C system headers
// includes for C++ system headers
// includes from other libraries
// includes from MCMC

namespace MarkovChainMonteCarlo
{

namespace Chain
{

/**
 * @class ChainBlock
 * @ingroup Chain
 * @brief Nodes of the linked list that represents the chain
 * @author James Till Matta
 * 
 * \tparam ParamType The floating point type to be used for the chain, float, double, long double, etc.
 * \tparam BlockSize The number of steps per walker that the block will hold
 * 
 * This class holds 'chunks' of the chain that is distributed across a linked
 * list style structure. Each block holds BlockSize parameter sets for all the
 * walkers
 * 
 */
template <class ParamType, int BlockSize>
class ChainBlock
{
public:
    /*!
     * \brief ChainBlock Constructor
     * \param prev Previous block in the chain, can be nullptr for the head node
     * \param numWalkers Number of walkers to store parameters for
     * \param numCellsPerWalker Number of parameters plus overhead that the walker needs
     */
    ChainBlock(ChainBlock<ParamType, BlockSize>* prev, int numWalkers, int numCellsPerWalker);
    ~ChainBlock(){delete[] chainArray;}
    
    /*!
     * \brief isFull checks if there are remaining chain steps available in the block
     * \return True if there are still steps open, false otherwise
     */
    bool isFull(){return (firstEmptyStep >= BlockSize);}
    
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
    bool incrementChainStepAndCheckNotFull(){++firstEmptyStep; return (firstEmptyStep < BlockSize);}
    
    friend class Chain;
    friend class ChainStepIterator;
    friend class ChainPsetIterator;
private:
    //Linked list book-keeping
    ChainBlock<ParamType, BlockSize>* lastBlock = nullptr; ///<pointer to the previous block in the chain linked list
    ChainBlock<ParamType, BlockSize>* nextBlock = nullptr; ///<pointer to the next block in the chain linked list
    //Chain book-keeping
    int walkerCount; ///<Number of walkers included in this chain
    int cellsPerWalker; ///<Number of cells needed by each walker
    int cellsPerStep; ///<Number of cells for all the walkers in a single step
    int firstEmptyStep = 0; ///<Index (in units of cellsPerStep) of the next empty step
    ParamType* chainArray = nullptr;///<Pointer into the array that stores the parameter set for every walker for every step in the block
};

template <class ParamType, int BlockSize>
ChainBlock<ParamType, BlockSize>::ChainBlock(ChainBlock<ParamType, BlockSize>* prev, int numWalkers, int numCellsPerWalker):
    lastBlock(prev), walkerCount(numWalkers), cellsPerWalker(numCellsPerWalker),
    cellsPerStep(numWalkers*numCellsPerWalker)
{
    chainArray = new ParamType[BlockSize*cellsPerStep];
}

template <class ParamType, int BlockSize>
void ChainBlock<ParamType, BlockSize>::storeWalker(int walkerNum, ParamType* walkerData)
{
    ParamType* offset = chainArray + (firstEmptyStep*cellsPerStep + walkerNum*cellsPerWalker);
    for(int i=0; i<cellsPerWalker; ++i)
    {
        offset[i] = walkerData[i];
    }
}

}
}
#endif  //MCMC_CHAIN_CHAINBLOCK_H
