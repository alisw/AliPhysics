////////////////////////////////////////////////////////////////////////////////
//
// Author: Artur Szostak
// Email:  artur@alice.phy.uct.ac.za | artursz@iafrica.com
//
////////////////////////////////////////////////////////////////////////////////

#ifndef dHLT_SYSTEM_SYSTEMERROR_HPP
#define dHLT_SYSTEM_SYSTEMERROR_HPP

#include "BasicTypes.hpp"
#include "Error.hpp"

namespace dHLT
{
namespace System
{


class Error : public dHLT::Error
{
public:

	Error() throw();

	Error(const Int code) throw()
	{
		errorcode = code; 
	};

	virtual const char* Message() const throw()
	{
		return AsString(errorcode);
	};
	
	virtual Int ErrorCode() const throw()
	{
		return errorcode;
	};
	
	static const char* AsString(const Int errorcode) throw();

protected:

	Int errorcode;
};


} // System
} // dHLT

#endif // dHLT_SYSTEM_SYSTEMERROR_HPP
