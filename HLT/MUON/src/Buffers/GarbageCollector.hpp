////////////////////////////////////////////////////////////////////////////////
//
// Author: Artur Szostak
// Email:  artur@alice.phy.uct.ac.za | artursz@iafrica.com
//
////////////////////////////////////////////////////////////////////////////////

#ifndef dHLT_BUFFERS_GARBAGE_COLLECTOR_HPP
#define dHLT_BUFFERS_GARBAGE_COLLECTOR_HPP

#include "BasicTypes.hpp"

namespace dHLT
{
namespace Buffers
{


class GarbageCollector
{
public:

	/* Allocates 'size' bytes of memory from the heap and sets the reference
	   count to 1.
	 */
	static void* Allocate(const UInt size);

	/* Increments the reference count. Every time the same block of memory is
	   referenced this method should be called.
	   WARNING! only use memory pointers returned by the Allocate method as
	   the parameter to IncRefCount.
	 */
	static void IncRefCount(void* memory);

	/* Decrements the reference count. If the count is zero then the memory
	   is released.
	 */
	static void Free(void* memory);
};


} // Buffers
} // dHLT

#endif // dHLT_BUFFERS_GARBAGE_COLLECTOR_HPP
