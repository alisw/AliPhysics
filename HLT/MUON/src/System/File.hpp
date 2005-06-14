////////////////////////////////////////////////////////////////////////////////
//
// Author: Artur Szostak
// Email:  artur@alice.phy.uct.ac.za | artursz@iafrica.com
//
////////////////////////////////////////////////////////////////////////////////

#ifndef dHLT_SYSTEM_FILE_HPP
#define dHLT_SYSTEM_FILE_HPP

#include "BasicTypes.hpp"
#include "System/SystemTypes.hpp"
#include "System/SystemError.hpp"

namespace dHLT
{
namespace System
{


class File
{
public:
	
	enum Mode
	{
		READ     = 0x01,
		WRITE    = 0x02,
		CREATE   = 0x04,
		TRUNCATE = 0x08,
		APPEND   = 0x10
	};
	
	enum Location
	{
		FromStart,
		FromEnd,
		FromCurrent
	};

	File() throw (System::Error);
	
	enum StdStream
	{
		StdIn,
		StdOut,
		StdError
	};

	File(const StdStream stream) throw (System::Error);
	
	/* Create a new file object and open the file.
	 */
	File(const char* filename, const UInt mode = READ | WRITE | CREATE) throw (System::Error);
	
	~File() throw (System::Error);
	
	/* Open the specified file. The mode must be bitwise-or'd with zero or more
	   of the following:
		READ      : File should be opened readable.
		WRITE     : File should be opened writeable.
		CREATE    : File should be created if it does not exist.
		TRUNCATE  : File should be truncated to zero.
		APPEND    : Set the file pointer to the end of the file.
	 */
	void Open(const char* filename, const UInt mode = READ | CREATE)
		throw (System::Error);

	/* Closes the file.
	 */
	void Close() throw (System::Error);
	
	/* Return the size of the open file.
	   The file must be open.
	 */
	UInt Size() const throw (System::Error);
	
	/* Checks if the open file is a regular file.
	 */
	bool IsARegularFile() const throw (System::Error);
	
	/* Moves the file pointer offset number of bytes relative to the specified location:
	      FromStart   : New file pointer = offset number of bytes.
	      FromEnd     : New file pointer = size of file + offset number of bytes.
	      FromCurrent : New file pointer = current file pointer + offset number of bytes.
	   The file must be open before calling this method.
	 */
	UInt Seek(UInt offset, const Location location = FromCurrent) throw (System::Error);
	
	/* Resize the file to the specified new size.
	   The file must be open before calling this method.
	 */
	void Resize(const UInt newsize) throw (System::Error);
	
	/* Truncate the file to the current Seek position.
	   The file must be open before calling this method.
	 */
	void Truncate() throw (System::Error);
	
	/* Reads 'size' number of bytes from the file and stores them into buffer.
	   The caller is responsible for allocating enough memory pointer to by buffer.
	   The actual number of bytes read is returned. This may be less than size.
	   The file must be open before calling this method.
	 */
	UInt Read(char* buffer, const UInt size) const throw (System::Error);
	UInt Read(void* buffer, const UInt size) const throw (System::Error);
	
	/* Writes 'size' number of bytes from buffer to the file. The actual number
	   of bytes written is returned. This may be less than size.
	   The file must be open before calling this method.
	 */
	UInt Write(const char* buffer, const UInt size) throw (System::Error);
	UInt Write(const void* buffer, const UInt size) throw (System::Error);
	
	// Return the system file handle.
	SystemFile Handle() const { return file; };

private:
	
	void OpenFile(const char* filename, const UInt mode);
	void CloseFile();

	SystemFile file;
};


} // System
} // dHLT

#endif // dHLT_SYSTEM_FILE_HPP
