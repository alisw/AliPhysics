////////////////////////////////////////////////////////////////////////////////
//
// Author: Artur Szostak
// Email:  artur@alice.phy.uct.ac.za | artursz@iafrica.com
//
////////////////////////////////////////////////////////////////////////////////

#ifndef dHLT_SYSTEM_DIRECTORY_HPP
#define dHLT_SYSTEM_DIRECTORY_HPP

#include "BasicTypes.hpp"
#include "System/SystemTypes.hpp"
#include "System/SystemError.hpp"

namespace dHLT
{
namespace System
{


class Directory
{
public:

	Directory() throw (System::Error);
	
	/* Create a directory class and open the specified directory for reading.
	 */
	Directory(const char* path) throw (System::Error);
	
	~Directory() throw (System::Error);
	
	/* Open the specified directory for reading.
	 */
	void Open(const char* path) throw (System::Error);

	/* Close the directory handle.
	 */
	void Close() throw (System::Error);
	
	/* Rewind the directory to the start. i.e. reset the state like it was just
	   after opening the directory.
	 */
	void Reset() const throw (System::Error);
	
	/* Read the next filename entry in the directory.
	   NOTE: This method is non-reentrant! Use the version below for multi-threaded
	   applications. Unless you can guarantee only one thread uses this method at
	   a time.
	 */
	const char* Read() const throw (System::Error);
	
	/* Reentrant version of the above read method. (only if getdents syscall is available)
	   Returns false if length of the filename buffer was not large
	   enough to contain the whole filename. The filename buffer is untouched
	   if false is returned.
	   An empty filename string indicates the end of directory listing.
	 */
	bool Read(char* filename, const UInt length) const throw (System::Error);
	
	SystemDirectory Handle() const { return dir; };

private:
	
	void OpenDir(const char* path);
	void CloseDir();

	SystemDirectory dir;
};


} // System
} // dHLT

#endif // dHLT_SYSTEM_DIRECTORY_HPP
