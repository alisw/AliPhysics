////////////////////////////////////////////////////////////////////////////////
//
// Author: Artur Szostak
// Email:  artur@alice.phy.uct.ac.za | artursz@iafrica.com
//
////////////////////////////////////////////////////////////////////////////////

#include "System/Routines.hpp"

#include <sys/types.h>
#include <sys/stat.h>
#include <unistd.h>
#include <asm/errno.h>

namespace dHLT
{
namespace System
{


UInt GetFileSize(const char* filename) throw (System::Error)
{
	struct stat buf;
	if ( stat(filename, &buf) != 0 ) throw System::Error();
	return buf.st_size;
};


bool IsARegularFile(const char* path) throw (System::Error)
{
	struct stat buf;
	int result = stat(path, &buf);
	if ( result == ENOENT ) return false;
	if ( result != 0 ) throw System::Error();
	return buf.st_mode & S_IFREG;
};


bool IsADirectory(const char* path) throw (System::Error)
{
	struct stat buf;
	int result = stat(path, &buf);
	if ( result == ENOENT ) return false;
	if ( result != 0 ) throw System::Error();
	return buf.st_mode & S_IFDIR;
};


} // System
} // dHLT
