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
#ifndef MCMC_UTILITY_MISC_H
#define MCMC_UTILITY_MISC_H
// includes for C system headers
// includes for C++ system headers
// includes from other libraries
// includes from MCMC

namespace MCMC
{

namespace Utility
{

static const unsigned int AlignmentLength = 64;///<Stores the memory boundary to force memory alignment to, 64 is sufficient for cache lines and up to the 256-bit AVX instructions, 128 will handle AVX-512 instructions as well

}
}
#endif  //MCMC_UTILITY_MISC_H
