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
#ifndef MCMC_UTILITY_ARRAYDELETER_H
#define MCMC_UTILITY_ARRAYDELETER_H
// includes for C system headers
// includes for C++ system headers
// includes from other libraries
// includes from MCMC

namespace MCMC
{

namespace Utility
{

/**
 * @class ArrayDeleter
 * @ingroup Utility
 * @brief A custom deleter for shared pointers that own arrays
 * @author James Till Matta
 * 
 * \tparam ArrayType The datatype the array is made of
 * 
 * This class simply provides the correct deleter to shared pointer. It is needed
 * for C++11 through C++14, in C++17 using the syntax std::shared_ptr<ParamType[]>
 * removes the need for this deleter
 */
template <class ArrayType>
class ArrayDeleter
{
public:
    void operator()(ArrayType const* arrayPtr){delete[] arrayPtr;}
private:
};

}
}
#endif  //MCMC_UTILITY_ARRAYDELETER_H
