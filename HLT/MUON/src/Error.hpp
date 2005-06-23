////////////////////////////////////////////////////////////////////////////////
//
// Author: Artur Szostak
// Email:  artur@alice.phy.uct.ac.za | artursz@iafrica.com
//
////////////////////////////////////////////////////////////////////////////////

#ifndef dHLT_ERROR_HPP
#define dHLT_ERROR_HPP

#include "BasicTypes.hpp"
#include <exception>
#include <Riostream.h>

namespace dHLT
{


class Error : public std::exception
{
public:

	Error() throw() {};
	virtual ~Error() throw() {};

	/* Should return a human readable string containing a description of the
	   error.
	 */
	virtual const char* Message() const throw() = 0;
	
	/* Returns an error code describing the error. The error code should be
	   unique to the entire system
	 */
	virtual Int ErrorCode() const throw() = 0;
	
	virtual const char* what() const throw()
	{
		return Message();
	};
	
	/* Define the << operator for streams to be able to do something like:

               Error myerror;
               cout << myerror << endl;
	*/
	friend ostream& operator << (ostream& os, const dHLT::Error& error)
	{
		os << error.Message();
		return os;
	};
};


class OutOfMemory : public Error
{
public:
	virtual const char* Message() const throw();
	virtual Int ErrorCode() const throw();
};


/* When one needs to indicate that no more memory is available one should use the
   ThrowOutOfMemory method rather than explicitly using the code
       throw OutOfMemory();
   This is because the ThrowOutOfMemory routine throws a preallocated object so
   we are safe from having to allocate more (nonexistant) memory.
 */
void ThrowOutOfMemory() throw (OutOfMemory);


// Error code declarations.
enum
{
	OUT_OF_MEMORY = 0x10000001
};


} // dHLT

#endif // dHLT_ERROR_HPP
