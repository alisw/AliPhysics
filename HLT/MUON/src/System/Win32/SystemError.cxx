////////////////////////////////////////////////////////////////////////////////
//
// Author: Artur Szostak
// Email:  artur@alice.phy.uct.ac.za | artursz@iafrica.com
//
////////////////////////////////////////////////////////////////////////////////

#include "System/SystemError.hpp"

#include <windows.h>

namespace dHLT
{
namespace System
{


Error::Error() throw() : dHLT::Error()
{
	errorcode = GetLastError();
};


const char* Error::AsString(const Int errorcode) throw()
{
	static char strbuf[1024];
	char* str = (char*) &strbuf[0];

	FormatMessage(
		FORMAT_MESSAGE_FROM_SYSTEM | FORMAT_MESSAGE_IGNORE_INSERTS,
		NULL,
		errorcode,
		MAKELANGID(LANG_NEUTRAL, SUBLANG_DEFAULT), // Default language
		(LPTSTR) str,
		sizeof(strbuf),
		NULL 
	);

	return str;
};


} // System
} // dHLT
