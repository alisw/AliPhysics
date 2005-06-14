////////////////////////////////////////////////////////////////////////////////
//
// Author: Artur Szostak
// Email:  artur@alice.phy.uct.ac.za | artursz@iafrica.com
//
////////////////////////////////////////////////////////////////////////////////

#include "Version.hpp"
#ifndef MAJOR_VERSION
 #include "VersionNumbers.hpp"
#endif
#include <stdio.h>

namespace dHLT
{

const Char* VersionString()
{
	static Char strbuf[40];
	Char* str = (Char*) &strbuf;
	sprintf(str, "%d.%d.%d", MAJOR_VERSION, MINOR_VERSION, BUILD_NUMBER);
	return str;
};

UInt MajorVersion()
{
	return MAJOR_VERSION;
};

UInt MinorVersion()
{
	return MINOR_VERSION;
};

UInt BuildNumber()
{
	return BUILD_NUMBER;
};


}; // dHLT
