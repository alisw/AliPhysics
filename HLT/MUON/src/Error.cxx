////////////////////////////////////////////////////////////////////////////////
//
// Author: Artur Szostak
// Email:  artur@alice.phy.uct.ac.za | artursz@iafrica.com
//
////////////////////////////////////////////////////////////////////////////////

/* AliHLTMUONCoreError is the base excpetion class used by the dHLT subsystem.
   All child classes used to throw exception should be derived from this
   class to allow easy catching of classes of errors.
   
   AliHLTMUONCoreOutOfMemory is also defined to be used when the system runs
   out of memory. Do not throw this object directly but rather use
   AliHLTMUONCoreThrowOutOfMemory which throws a pree allocated static object.
 */
 
#include "Error.hpp"
#include "Utils.hpp"

namespace
{
	// The one and only pree allocated out of memory error object.
	static AliHLTMUONCoreOutOfMemory gAliOutOfMemObject;

} // end of namespace


const char* AliHLTMUONCoreOutOfMemory::Message() const throw()
{
	return "Out of memory.";
}

Int AliHLTMUONCoreOutOfMemory::ErrorCode() const throw()
{
	return kOutOfMemory;
}

void AliHLTMUONCoreThrowOutOfMemory() throw (AliHLTMUONCoreOutOfMemory)
{
	throw gAliOutOfMemObject;
}

