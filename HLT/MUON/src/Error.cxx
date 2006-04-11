////////////////////////////////////////////////////////////////////////////////
//
// Author: Artur Szostak
// Email:  artur@alice.phy.uct.ac.za | artursz@iafrica.com
//
////////////////////////////////////////////////////////////////////////////////

#include "Error.hpp"
#include "Utils.hpp"

namespace
{
	// The one and only pree allocated out of memory error object.
	static AliHLTMUONCoreOutOfMemory _out_of_memory_;

} // end of namespace


const char* AliHLTMUONCoreOutOfMemory::Message() const throw()
{
	return "Out of memory.";
}

Int AliHLTMUONCoreOutOfMemory::ErrorCode() const throw()
{
	return kOUT_OF_MEMORY;
}

void AliHLTMUONCoreThrowOutOfMemory() throw (AliHLTMUONCoreOutOfMemory)
{
	throw _out_of_memory_;
}

