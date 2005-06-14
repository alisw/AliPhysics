////////////////////////////////////////////////////////////////////////////////
//
// Author: Artur Szostak
// Email:  artur@alice.phy.uct.ac.za | artursz@iafrica.com
//
////////////////////////////////////////////////////////////////////////////////

#include "System/SystemError.hpp"

#include <errno.h>
#include <string.h>

namespace dHLT
{
namespace System
{


Error::Error() throw() : dHLT::Error()
{
	errorcode = errno;
};

const char* Error::AsString(const Int errorcode) throw()
{
	return strerror(errorcode);
};


} // System
} // dHLT
