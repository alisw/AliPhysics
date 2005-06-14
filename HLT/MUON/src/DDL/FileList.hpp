////////////////////////////////////////////////////////////////////////////////
//
// Author: Artur Szostak
// Email:  artur@alice.phy.uct.ac.za | artursz@iafrica.com
//
// Author: Gareth de Vaux
// Email:  devaux@lhc.phy.uct.ac.za | dhlt@lordcow.org
//
////////////////////////////////////////////////////////////////////////////////

#ifndef dHLT_DDL_FILE_LIST_HPP
#define dHLT_DDL_FILE_LIST_HPP

#include "Error.hpp"
#include "BasicTypes.hpp"

#include <vector>

namespace dHLT
{
namespace DDL
{


enum
{
	PATH_NOT_FOUND = 0x10020001
};


class PathNotFound : public Error
{
public:
	PathNotFound(const char* path) throw (OutOfMemory);
	virtual ~PathNotFound() throw ();

	virtual const char* Message() const throw();
	virtual Int ErrorCode() const throw();
	
private:

	char* message;
};


class FileList
{
public:

	FileList();
	
	/* Creates a new file list and initialises it by adding the specified path.
	 */
	FileList(const char* path, const bool recursive = true);
	
	~FileList();
	
	/* Returns true if the file list is empty.
	 */
	bool Empty() const;
	
	/* Returns the number of files in the list.
	 */
	UInt Count() const;

	/* If path contains a valid directory then all files in that directory are
	   added to the filelist. If recursive is true then all files in the sub
	   directories are also added.
	   If the path contains the name of a specific file then it is assumed the
	   file is a text file containing a list of paths to add to the list.
	   In this case the recursive flag should be applied to all the paths found.
	 */ 
	void Add(const char* path, const bool recursive = true);
	
	/* Checks if the file list contains the specified file name.
	 */
	bool Contains(const char* filename) const;
	
	/* Looks for the specifed file name in the file list and returns the index
	   number to that file. 
	   -1 is returned if the file could not be found.
	 */
	Int Find(const char* filename) const;
	
	/* Returns the index'th file name in the list.
	 */
	const char* FileName(const UInt index) const;
	
	/* Returns the size of the index'th file in the list.
	 */
	UInt FileSize(const UInt index) const;
	
	/* Loads the specified file contents into memory. The calling process must 
	   clean up the memory allocated with a call to delete [] buffer.
	   size is filled with the size of the buffer.
	   True is returned if the file was located in the file list.
	 */
	bool Fetch(const char* filename, char*& buffer, UInt& size) const;
	
	// Same as above but fetch the index'th file.
	void Fetch(const UInt index, char*& buffer, UInt& size) const;
	
	/* Loads the file contents just like above, however the buffer must be
	   preallocated by the caller. Only 'size' number of bytes are read from
	   the file and filled into the buffer. If the file is smaller than the
	   buffer size then the remaining bytes are left untouched.
	 */
	bool Fetch(const char* filename, const UInt size, char* buffer) const;
	
	// Same as above but fetch the index'th file.
	void Fetch(const UInt index, const UInt size, char* buffer) const;
	
	/* Empties the file list.
	 */
	void Clear();

private:

	void AddFile(const char* filename);  // Called by Add & constructor

	void AddDir(const char* path, const bool recursive = true);

	UInt ParseLine(const UInt currentStart, char* buffer, const UInt end); 

	std::vector<char*> list;
};



} // DDL
} // dHLT

#endif // dHLT_DDL_FILE_LIST_HPP
