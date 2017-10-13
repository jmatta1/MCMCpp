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
#ifndef MCMC_UTILITY_USEROBJECTTEST_H
#define MCMC_UTILITY_USEROBJECTTEST_H
// includes for C system headers
// includes for C++ system headers
#include<type_traits>
// includes from other libraries
// includes from MCMC

namespace MarkovChainMonteCarlo
{

namespace Utility
{

namespace Detail
{
//Check Type Equality Object
//default case, the two types are not equal
/*!
 * @struct SameType
 * @brief This is a simple class using template specialization to figure out
 * if two classes are the same time or not
 * 
 * @tparam T1 The first object type to be compared
 * @tparam T2 The second ojbect type to be compared
 */
template<class T1, class T2>
struct SameType { typedef std::false_type ValType; };

/*!
 * @cond HIDDEN_SYMBOLS
 */
template<class T1>
struct SameType<T1, T1> { typedef std::true_type ValType;};
/*!
 * @endcond
 */


/*!
 * @brief testSignature_PerformAction tests the existance of a function named PerformAction that has 1 argument
 * @tparam TestClass The class to test for existence of a member function
 * @tparam RetType The type of the return argument
 * @tparam Arg0 The type of the first function argument
 */
template<class TestClass, class RetType, class Arg0>
static auto testSignature_performAction(unsigned long long) -> typename SameType<RetType, decltype(std::declval<TestClass>().performAction(std::declval<Arg0>())) >::ValType;
/*!
 * @cond HIDDEN_SYMBOLS
 * substitution failure branch, either the function does not exist or the arguments passed cannot be coerced into being the correct arguments
 */
template<class, class, class >
static auto testSignature_performAction(int) -> std::false_type;
/*!
 * @endcond
 */

/*!
 * @brief testSignature_PerformAction tests the existance of operator() that has 1 argument
 * @tparam TestClass The class to test for existence of a member function
 * @tparam RetType The type of the return argument
 * @tparam Arg0 The type of the first function argument
 */
template<class TestClass, class RetType, class Arg0>
static auto testSignature_functor(unsigned long long) -> typename SameType<RetType,
    decltype(std::declval<TestClass>().operator()(std::declval<Arg0>()))>::ValType;
/*!
 * @cond HIDDEN_SYMBOLS
 * substitution failure branch, either the function does not exist or the arguments passed cannot be coerced into being the correct arguments
 */
template<class, class, class >
static auto testSignature_functor(int) -> std::false_type;
/*!
 * @endcond
 */
}

/*!
 * @brief Outer Function to perform test for perform action function
 */
template<class TestClass, class RetVal, class Arg0>
struct CheckPerformAction : decltype(Detail::testSignature_performAction<TestClass, RetVal, Arg0>(0ULL) ){};

/*!
 * @brief Outer Function to perform test for operator() function
 */
template<class TestClass, class RetVal, class Arg0>
struct CheckFunctor : decltype( Detail::testSignature_functor<TestClass, RetVal, Arg0>(0ULL) ){};

}
}
#endif  //MCMC_UTILITY_USEROBJECTTEST_H
