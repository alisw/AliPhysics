////////////////////////////////////////////////////////////////////////////////
//
// Author: Artur Szostak
// Email:  artur@alice.phy.uct.ac.za | artursz@iafrica.com
//
////////////////////////////////////////////////////////////////////////////////

#include "Buffers/GarbageCollector.hpp"
#include <stdlib.h>
#include "Utils.hpp"

#ifdef DEBUG

#include <vector>
#include <algorithm>

namespace
{
	std::vector<void*> pointerlist;
	
	/* Adds a new pointer to the pointer list.
	   This routine is called from GarbageCollector::Allocate when a new
	   memory block is allocated to allow us to check if all memory was
	   cleaned up properly.
	 */
	void AddNewPointer(void* p)
	{
		pointerlist.push_back(p);
	};
	
	/* Removes the pointer p from the pointer list. If it is not found then
	   an assertion failure is generated.
	 */
	void RemovePointer(void* p)
	{
		std::vector<void*>::iterator elem = 
			find(pointerlist.begin(), pointerlist.end(), p);
		Assert( elem != pointerlist.end() );
		pointerlist.erase(elem);
	};
	
	
#	ifndef NO_ON_CLEANUP_CHECK

	/* This OnCleanup class is used to check if the pointer list is empty
	   at the end of a program. (i.e. that all memory was released properly)
	   Note: since some destructors for global objects might only be called
	   after this objects destructor, one might need to reorder the object
	   files during linking to eliminate false checks on correct code.
	   If that does not work and global instanciations of object can not be
	   moved to a more local scope such as the main routine then one should
	   just skip this test or make it explicit, by defining NO_ON_CLEANUP_CHECK
	 */
	class OnCleanup
	{
	public:
		~OnCleanup()
		{
			// Check to see if all memory was released.
			Assert( pointerlist.size() == 0 );
		};
	} oncleanup;

#	endif // NO_ON_CLEANUP_CHECK
	
};

#endif // DEBUG


namespace dHLT
{
namespace Buffers
{


// The memory layout of allocated blocks is as follows:
// 
//  Start of memory block.
//  | 
//  V
//  [ Ref count ] [ User data .... ]
//                ^
//                |
//               pointer returned to user.
//
// The reference count is an unsigned integer located just before the start of
// the users memory region. Only if this ref count reaches zero will the memory
// block be deallocated.


void* GarbageCollector::Allocate(const UInt size)
{
	UInt* block = (UInt*) malloc(size + sizeof(UInt));
	void* usermemory = (void*)(block + 1);
	*block = 1;
	DebugCode( AddNewPointer(usermemory); );
	return usermemory;
};


void GarbageCollector::IncRefCount(void* memory)
{
	UInt* block = ((UInt*)memory) - 1;
	(*block)++;
};


void GarbageCollector::Free(void* memory)
{
	UInt* block = ((UInt*)memory) - 1;
	(*block)--;
	if (*block == 0)
	{
		DebugCode( RemovePointer(memory); );
		free(block);
	};
};


} // Buffers
} // dHLT
