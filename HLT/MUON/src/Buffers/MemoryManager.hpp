////////////////////////////////////////////////////////////////////////////////
//
// Author: Artur Szostak
// Email:  artur@alice.phy.uct.ac.za | artursz@iafrica.com
//
////////////////////////////////////////////////////////////////////////////////

#ifndef dHLT_BUFFERS_MEMORY_MANAGER_HPP
#define dHLT_BUFFERS_MEMORY_MANAGER_HPP

#include "BasicTypes.hpp"

namespace dHLT
{
namespace Buffers
{


template <typename DataType, UInt size = 1024>
class MemoryManager
{
public:

	MemoryManager()
	{
		current = &data[0];
	};


	~MemoryManager()
	{
	};


	void* Allocate()
	{
		return current++;
		//return malloc( sizeof(DataType));
	};


	void Free(void* memory)
	{
		//free(memory);
	};


	void Compress()
	{
		
	};


private:
	
	struct Block
	{
		DataType data[size];
	};

	DataType* current;
	DataType data[1000000];

};


} // Buffers
} // dHLT

#endif // dHLT_BUFFERS_MEMORY_MANAGER_HPP
