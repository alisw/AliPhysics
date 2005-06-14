////////////////////////////////////////////////////////////////////////////////
//
// Author: Artur Szostak
// Email:  artur@alice.phy.uct.ac.za | artursz@iafrica.com
//
////////////////////////////////////////////////////////////////////////////////

#ifndef dHLT_BASIC_TYPES_HPP
#define dHLT_BASIC_TYPES_HPP

namespace dHLT
{


#ifndef NULL
#	define NULL 0x0
#endif


typedef char Char;
typedef short Short;
typedef int Int;
typedef long Long;
typedef unsigned char UChar;
typedef unsigned short UShort;
typedef unsigned int UInt;
typedef unsigned long ULong;
typedef float Float;
typedef double Double;

typedef char* NullString;


} // dHLT

#endif // dHLT_BASIC_TYPES_HPP
