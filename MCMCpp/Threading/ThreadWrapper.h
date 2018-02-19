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
#ifndef MCMCPP_THREADING_THREADWRAPPER_H
#define MCMCPP_THREADING_THREADWRAPPER_H
// includes for C system headers
// includes for C++ system headers
// includes from other libraries
// includes from MCMCpp

namespace MCMC
{
namespace Threading
{
/*!
 * @class ThreadWrapper
 * @ingroup Threading
 * @brief Simple class allowing complex thread objects to be wrapped and used without copying
 * @author James Till Matta
 * 
 * @tparam ThreadType The class representing the thread object
 * 
 * Because std::thread copies the callable object that serves as the thread of
 * execution can cause problems for large and complex objects it is better to
 * make a simple wrapper that can be used to wrap a pointer to the object and
 * then that is passed to the thread object. This means that all that gets
 * copied is a pointer minimizing insanity and redundant copies of objects
 */
template<class ThreadType>
struct ThreadWrapper
{
public:
    /*!
     * \brief ThreadWrapper The standard constructor of this thread wrapper
     * \param p A pointer to the object that contains the actual thread of execution
     */
    ThreadWrapper(ThreadType* p): ptr(p){}
    
    /*!
     * \brief ThreadWrapper The copy constructor, generates a copy of this class
     * \param copy The ThreadWrapper class to be copied
     */
    ThreadWrapper(const ThreadWrapper<ThreadType>& copy): ptr(copy.ptr){}
    /*!
     * @brief ~ThreadWrapper Allow ptr to destruct without deleting its contents, something else owns the contents of ptr
     */
    ~ThreadWrapper(){}
    
    /*!
     * \brief operator= An assignement operator to go with the copy constructor, doesn't delete original contents, something else owns that
     * \param rhs The ThreadWrapper being assigned to this object
     * \return A reference to this object
     */
    ThreadWrapper<ThreadType>& operator=(const ThreadWrapper<ThreadType>& rhs){ptr = rhs.ptr; return *this;}
    
    /*!
     * \brief operator() Allows the thread object to call the underlying ThreadType object by forwarding the call
     */
    void operator()() { (*ptr)(); }
    
private:
    ThreadType* ptr; ///<A pointer to the actual thread object
};

}
}
#endif  //MCMCPP_THREADING_THREADWRAPPER_H
