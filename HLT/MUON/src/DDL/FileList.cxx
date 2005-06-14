///////////////////////////////////////////////////////////////////////////////
//
// Author: Artur Szostak
// Email:  artur@alice.phy.uct.ac.za | artursz@iafrica.com
//
// Author: Gareth de Vaux
// Email:  devaux@lhc.phy.uct.ac.za | dhlt@lordcow.org
//
////////////////////////////////////////////////////////////////////////////////

#include "FileList.hpp"
#include "Utils.hpp"
#include "System/File.hpp"
#include "System/Directory.hpp"
#include "System/Routines.hpp"

#include <stdlib.h>
#include <string.h>
#include <string>

using std::string;

namespace dHLT
{
namespace DDL
{


PathNotFound::PathNotFound(const char* path) throw (OutOfMemory)
{
	char* str = "Could not find path: ";
	UInt strsize = strlen(str);
	UInt pathsize = strlen(str);
	message = new char[strsize + pathsize + 1];
	memcpy(message, str, strsize);
	memcpy(&message[strsize], path, pathsize);
	message[strsize + pathsize] = '\0';  // Set the NULL character at end of string.
};


PathNotFound::~PathNotFound() throw ()
{
	Assert( message != NULL );
	delete [] message;
};


const char* PathNotFound::Message() const throw()
{
	return message;
};


Int PathNotFound::ErrorCode() const throw()
{
	return PATH_NOT_FOUND;
};

////////////////////////////////////////////////////////////////////////////////

using System::File;
using System::Directory;


FileList::FileList()
{
	// Nothing special to be done.
};

FileList::FileList(const char* path, const bool recursive)
{
	Add(path, recursive);
};

FileList::~FileList()
{
	Clear();
};

bool FileList::Empty() const
{
	return list.empty();
};

UInt FileList::Count() const
{
	return list.size();
};

void FileList::Add(const char* path, const bool recursive)
{
	if (System::IsADirectory(path))
		AddDir(path, recursive);
	else
	{
		if (System::IsARegularFile(path))
		{
			UInt size, bytesread = 0, index = 0;

			size = System::GetFileSize(path);
			char* buffer = new char[size];
			try
			{
				File file(path, File::READ);

				while (bytesread < size)
					bytesread += file.Read(buffer + bytesread, size - bytesread);

				while (index < size)
					index = ParseLine(index, buffer, size);
			}
			finally
			(
				delete [] buffer;
			);
		}
		else
			throw PathNotFound(path);
	};
};


bool FileList::Contains(const char* filename) const
{
	for (UInt i = 0; i < Count(); i++)
	{
		if ( strcmp(list[i], filename) == 0 )
			return true;

	};
	return false;
};


Int FileList::Find(const char* filename) const
{
	for (UInt i = 0; i < Count(); i++)
	{
		if ( strcmp(list[i], filename) == 0 )
		{
			Assert( (UInt) ((Int) i) == i );  // Check typecast conversion.
			return (Int) i;
		};
	};
	return -1;
};


const char* FileList::FileName(const UInt index) const
{
	Assert( index < Count() );
	return list[index];
};


UInt FileList::FileSize(const UInt index) const
{
	Assert( index < Count() );
	return System::GetFileSize( FileName(index) );
};


bool FileList::Fetch(const char* filename, char*& buffer, UInt& size) const
{
	Int index = Find(filename);
	if (index != -1)
	{
		Fetch(index, buffer, size);
		return true;
	}
	else
		return false;
};


void FileList::Fetch(const UInt index, char*& buffer, UInt& size) const
{
	Assert( index < Count() );
	
	const char* filename = FileName(index);
	size = System::GetFileSize(filename);
	buffer = new char[size];
	
	File file(filename, File::READ);
	UInt bytesread = 0;
	while (bytesread < size)
	{
		bytesread += file.Read(buffer + bytesread, size - bytesread);
	};
};


bool FileList::Fetch(const char* filename, const UInt size, char* buffer) const
{
	Int index = Find(filename);
	if (index != -1)
	{
		Fetch(index, size, buffer);
		return true;
	}
	else
		return false;
};


void FileList::Fetch(const UInt index, const UInt size, char* buffer) const
{
	Assert( index < Count() );
	
	const char* filename = FileName(index);
	UInt filesize = System::GetFileSize(filename);
	
	// Make sure we do not overflow the buffer.
	if (size < filesize)
		filesize = size;
	
	File file(filename, File::READ);
	UInt bytesread = 0;
	while (bytesread < filesize)
	{
		bytesread += file.Read(buffer + bytesread, filesize - bytesread);
	};
};


void FileList::Clear()
{
	// Release memory allocated to filenames then clear the vector array itself.
	for (UInt i = 0; i < Count(); i++)
	{
		char* filename = list[i];
		delete [] filename;
	};
	list.clear();
};


void FileList::AddFile(const char* filename)
{
	char* name = new char[strlen(filename) + 1];
	strcpy(name, filename);
	list.push_back(name);
};


void FileList::AddDir(const char* path, const bool recursive)
{
	const char* tmpFile;
	string stringPath(path), slash("/"), stringFile, tmpPath;         

	Directory dir = path;

	tmpFile = dir.Read(); // get rid've the . & .. files
	tmpFile = dir.Read();

	for (tmpFile = dir.Read(); tmpFile != NULL; tmpFile = dir.Read())
	{
		stringFile = tmpFile;
		tmpPath = stringPath + slash + stringFile;

		if (System::IsADirectory(tmpPath.c_str()) && recursive)
			Add(tmpPath.c_str(), recursive);
		else
		{
			if (System::IsARegularFile(tmpPath.c_str()))
				AddFile(tmpPath.c_str());
			else
				throw PathNotFound(tmpPath.c_str());
		}
	};

	dir.Close();
};


UInt FileList::ParseLine(const UInt currentStart, char* buffer, const UInt end)
{
	UInt length = currentStart, start = currentStart, finish, i, j;
	char *tmpLine;

	// Find a while line.
	while (buffer[length] != '\n' and length < end)
		length++;

	finish = length;

	// Remove whitespace from the beginning and end of the line.
	while ((buffer[start] == ' ' or buffer[start] == '\t' or buffer[start] == '\r') and start < length)
		start++;

	while ((buffer[finish] == ' ' or buffer[finish] == '\t' or buffer[finish] == '\r') and finish > currentStart)
		finish--;

	length = finish - start + 1;
	tmpLine = new char[length];
	try
	{
		for (i = start, j = 0; i < finish ; i++, j++)
			tmpLine[j] = buffer[i];

		tmpLine[j] = '\0';

		if (System::IsADirectory(tmpLine))
			AddDir(tmpLine);
		else if (System::IsARegularFile(tmpLine))
			AddFile(tmpLine);
		else
			throw PathNotFound(tmpLine);
	}
	finally
	(
		delete [] tmpLine;
	);
	return finish+1;
};


} // DDL
} // dHLT
