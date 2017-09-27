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

namespace MarkovChainMonteCarlo
{

namespace Chain
{

/**
 * @class Chain
 * @ingroup Chain
 * @brief The main access point of the chain which manages the linked list
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
    /*!
     * \brief Chain Constructor
     * \param numWalkers Number of walkers to store parameters for
     * \param numCellsPerWalker Number of parameters plus overhead that the walker needs
     * \param maxSize Maximum size in bytes that the chain is allowed to occupy (does not include overhead of nodes etc)
     */
    Chain(int numWalkers, int numCellsPerWalker, unsigned long long maxSize);
    ~Chain();
    
private:
    //Linked list book-keeping
    ChainBlock<ParamType, BlockSize>* head = nullptr; ///<pointer to the previous block in the chain linked list
    ChainBlock<ParamType, BlockSize>* curr = nullptr; ///<pointer to the next block in the chain linked list
    //Chain book-keeping
    int walkerCount; ///<Number of walkers included in this chain
    int cellsPerWalker; ///<Number of cells needed by each walker
    int maxBlocks; ///<Calculated maximum number of blocks to allocate
    int firstEmptyStep = 0; ///<Index (in units of cellsPerStep) of the next empty step
};


}
}
#endif  //MCMC_CHAIN_CHAIN_H
