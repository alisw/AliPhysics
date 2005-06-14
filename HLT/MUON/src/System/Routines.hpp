////////////////////////////////////////////////////////////////////////////////
//
// Author: Artur Szostak
// Email:  artur@alice.phy.uct.ac.za | artursz@iafrica.com
//
////////////////////////////////////////////////////////////////////////////////

#ifndef dHLT_SYSTEM_ROUTINES_HPP
#define dHLT_SYSTEM_ROUTINES_HPP

#include "BasicTypes.hpp"
#include "System/SystemError.hpp"

namespace dHLT
{
namespace System
{


/* Returns the size of the specified file in bytes.
 */
UInt GetFileSize(const char* filename) throw (System::Error);

/* Returns true if the path is regular file.
 */
bool IsARegularFile(const char* path) throw (System::Error);

/* Returns true if the path is a directory name.
 */
bool IsADirectory(const char* path) throw (System::Error);


} // System
} // dHLT

#endif // dHLT_SYSTEM_ROUTINES_HPP
