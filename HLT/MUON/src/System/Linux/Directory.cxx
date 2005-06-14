////////////////////////////////////////////////////////////////////////////////
//
// Author: Artur Szostak
// Email:  artur@alice.phy.uct.ac.za | artursz@iafrica.com
//
////////////////////////////////////////////////////////////////////////////////

#include "System/Directory.hpp"

#include <sys/types.h>
#include <sys/stat.h>
#include <unistd.h>
#include <fcntl.h>
#include <errno.h>

#ifdef USE_GETDENTS_SYSCALL
#	include <linux/types.h>
#	include <linux/dirent.h>
#	include <linux/unistd.h>

// If we dont have getdents defined then just use the syscall directly.
#	ifndef getdents
#		include <sys/syscall.h>
#		define getdents(fd, entry, size) syscall(SYS_getdents, fd, entry, size)
#	endif // getdents

#else // USE_GETDENTS_SYSCALL
#	include <dirent.h>
#endif // USE_GETDENTS_SYSCALL

namespace dHLT
{
namespace System
{


Directory::Directory() throw (System::Error)
{
#ifdef USE_GETDENTS_SYSCALL
	dir = -1;
#else // USE_GETDENTS_SYSCALL
	dir = NULL;
#endif // USE_GETDENTS_SYSCALL
};


Directory::Directory(const char* path) throw (System::Error)
{
	OpenDir(path);
};


Directory::~Directory() throw (System::Error)
{
#ifdef USE_GETDENTS_SYSCALL
	if (dir != -1) CloseDir();
#else // USE_GETDENTS_SYSCALL
	if (dir != NULL) CloseDir();
#endif // USE_GETDENTS_SYSCALL
};


void Directory::Open(const char* path) throw (System::Error)
{
#ifdef USE_GETDENTS_SYSCALL
	if (dir != -1) CloseDir();
#else // USE_GETDENTS_SYSCALL
	if (dir != NULL) CloseDir();
#endif // USE_GETDENTS_SYSCALL
	OpenDir(path);
};


void Directory::Close() throw (System::Error)
{
#ifdef USE_GETDENTS_SYSCALL
	if (dir != -1)
	{
		CloseDir();
		dir = -1;
	};
#else // USE_GETDENTS_SYSCALL
	if (dir != NULL)
	{
		CloseDir();
		dir = NULL;
	};
#endif // USE_GETDENTS_SYSCALL
};


void Directory::Reset() const throw (System::Error)
{
#ifdef USE_GETDENTS_SYSCALL
	off_t result = lseek(dir, 0, SEEK_SET);
	if (result == -1) throw System::Error();
#else // USE_GETDENTS_SYSCALL
	rewinddir(dir);
#endif // USE_GETDENTS_SYSCALL
};


const char* Directory::Read() const throw (System::Error)
{
#ifdef USE_GETDENTS_SYSCALL

	static struct dirent entry;
	int result = getdents(dir, &entry, sizeof(entry));
	if (result == -1) throw System::Error();
	if (result == 0) return NULL;
	
	// Move to the next directory entry.
	off_t newpos = lseek(dir, entry.d_off, SEEK_SET);
	if (newpos == -1) throw System::Error();
	return &entry.d_name[0];
	
#else // USE_GETDENTS_SYSCALL

	struct dirent* entry = readdir(dir);
	if (entry == NULL)
	{
		if (errno != 0)
			throw System::Error();
		else
			return NULL;
	}
	else
		return &entry->d_name[0];
	
#endif // USE_GETDENTS_SYSCALL
};


bool Directory::Read(char* filename, const UInt length) const throw (System::Error)
{
#ifdef USE_GETDENTS_SYSCALL

	struct dirent entry;
	int result = getdents(dir, &entry, sizeof(entry));
	if (result == -1) throw System::Error();
	if (result == 0)
	{
		if (length >= 1)
		{
			filename[0] = '\0';
			return true;
		}
		else
			return false;
	};
	
	// Move to the next directory entry.
	off_t newpos = lseek(dir, entry.d_off, SEEK_SET);
	if (newpos == -1) throw System::Error();
	
	char* str = &entry.d_name[0];
		
#else // USE_GETDENTS_SYSCALL

	struct dirent* entry = readdir(dir);
	if (entry == NULL)
	{
		if (errno != 0) throw System::Error();
		if (length >= 1)
		{
			filename[0] = '\0';
			return true;
		}
		else
			return false;
	};

	char* str = &entry->d_name[0];
	
#endif // USE_GETDENTS_SYSCALL

	// Fill the filename buffer and return.
	if( strlen(str) + 1 <= length )
	{
		strcpy(filename, str);
		return true;
	}
	else
		return false;
};


void Directory::OpenDir(const char* path)
{
#ifdef USE_GETDENTS_SYSCALL
	dir = open(path, O_RDONLY | O_DIRECTORY);
	if (dir == -1) throw System::Error();
#else // USE_GETDENTS_SYSCALL
	dir = opendir(path);
	if (dir == NULL) throw System::Error();
#endif // USE_GETDENTS_SYSCALL
};


void Directory::CloseDir()
{
#ifdef USE_GETDENTS_SYSCALL
	if ( close(dir) != 0 ) throw System::Error();
#else // USE_GETDENTS_SYSCALL
	if ( closedir(dir) != 0 ) throw System::Error();
#endif // USE_GETDENTS_SYSCALL
};


} // System
} // dHLT
