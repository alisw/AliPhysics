////////////////////////////////////////////////////////////////////////////////
//
// Author: Artur Szostak
// Email:  artur@alice.phy.uct.ac.za | artursz@iafrica.com
//
////////////////////////////////////////////////////////////////////////////////

#include "System/File.hpp"

#include <sys/types.h>
#include <sys/stat.h>
#include <unistd.h>
#include <fcntl.h>
#include <asm/errno.h>

namespace dHLT
{
namespace System
{


File::File() throw (System::Error)
{
	file = -1;
};


File::File(const StdStream stream) throw (System::Error)
{
	switch (stream)
	{
	case StdIn:    file = 0; break;
	case StdOut:   file = 1; break;
	case StdError: file = 2; break;
	default:       throw System::Error(EBADF);
	};
};


File::File(const char* filename, const UInt mode) throw (System::Error)
{
	OpenFile(filename, mode);
};


File::~File() throw (System::Error)
{
	if (file != -1)
		CloseFile();
};


void File::Open(const char* filename, const UInt mode) throw (System::Error)
{
	if (file != -1) CloseFile();
	OpenFile(filename, mode);
};


void File::Close() throw (System::Error)
{
	if (file != -1)
	{
		CloseFile();
		file = -1;
	};
};


UInt File::Size() const throw (System::Error)
{
	struct stat buf;
	if ( fstat(file, &buf) != 0 ) throw System::Error();
	return buf.st_size;
};


bool File::IsARegularFile() const throw (System::Error)
{
	struct stat buf;
	if ( fstat(file, &buf) != 0 ) throw System::Error();
	return buf.st_mode & S_IFREG;
};


UInt File::Seek(UInt offset, const Location location) throw (System::Error)
{
	int whence = -1;
	switch (location)
	{
	case FromStart:   whence = SEEK_SET; break;
	case FromEnd:     whence = SEEK_END; break;
	case FromCurrent: whence = SEEK_CUR; break;
	};
	off_t result = lseek(file, offset, whence);
	if (result == -1) throw System::Error();
	return result;
};


void File::Resize(const UInt newsize) throw (System::Error)
{
	if (ftruncate(file, newsize) != 0) throw System::Error();
};


void File::Truncate() throw (System::Error)
{
	UInt newsize = Seek(0);
	Resize(newsize);
};


UInt File::Read(char* buffer, const UInt size) const throw (System::Error)
{
	ssize_t result = read(file, (void*)buffer, size);
	if (result == -1) throw System::Error();
	return result;
};


UInt File::Read(void* buffer, const UInt size) const throw (System::Error)
{
	return Read((char*)buffer, size);
};


UInt File::Write(const char* buffer, const UInt size) throw (System::Error)
{
	ssize_t result = write(file, (void*)buffer, size);
	if (result == -1) throw System::Error();
	return result;
};


UInt File::Write(const void* buffer, const UInt size) throw (System::Error)
{
	return Write((char*)buffer, size);
};


void File::OpenFile(const char* filename, const UInt mode)
{
	int flags = 0;
	if (mode & READ)
	{
		if (mode & WRITE)
			flags |= O_RDWR;
		else
			flags |= O_RDONLY;
	}
	else
	{
		if (mode & WRITE) flags |= O_WRONLY;
	};
	
	if (mode & CREATE) flags |= O_CREAT;
	if (mode & TRUNCATE) flags |= O_TRUNC;
	
	file = open(filename, flags, S_IRUSR | S_IWUSR | S_IRGRP | S_IWGRP | S_IROTH);
	if (file == -1) throw System::Error();
	
	if (mode & APPEND)
		Seek(0, FromEnd);
};


void File::CloseFile()
{
	if ( 0 <= file and file <= 2 ) return;  // Do not touch the standard stream file handles.
	if ( close(file) != 0 ) throw System::Error();
};


} // System
} // dHLT
