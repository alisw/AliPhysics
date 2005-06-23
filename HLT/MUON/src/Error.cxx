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
	static dHLT::OutOfMemory _out_of_memory_;

} // end of namespace


namespace dHLT
{


const char* OutOfMemory::Message() const throw()
{
	return "Out of memory.";
}

Int OutOfMemory::ErrorCode() const throw()
{
	return OUT_OF_MEMORY;
}

void ThrowOutOfMemory() throw (OutOfMemory)
{
	throw _out_of_memory_;
}


} // dHLT
