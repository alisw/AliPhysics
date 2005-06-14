////////////////////////////////////////////////////////////////////////////////
//
// Author: Artur Szostak
// Email:  artur@alice.phy.uct.ac.za | artursz@iafrica.com
//
////////////////////////////////////////////////////////////////////////////////

#ifndef dHLT_BUFFERS_QUEUE_HPP
#define dHLT_BUFFERS_QUEUE_HPP

#include "BasicTypes.hpp"
#include "new.hpp"

namespace dHLT
{
namespace Buffers
{


template <typename DataType, UInt blocksize = 1024>
class Queue
{
public:

	Queue()
	{
		firstblock = lastblock = new Block( new Block() );
		firstblock->next->next = firstblock;
		firstblockend = lastblockend = (void*)( (UInt)firstblock + sizeof(Block) );
		firstelement = lastelement = &firstblock->data[0];
	};


	~Queue()
	{
		Block* current = firstblock->next;
		firstblock->next = NULL;  // break the cyclic block structure.
		Block* blocktodelete;
		while (current != NULL)
		{
			blocktodelete = current;
			current = current->next;
			delete blocktodelete;
		};
	};


	void Push(const DataType& data)
	{
		*lastelement = data;
		lastelement++;

		// If we hit the end of the block then we need to move the pointers
		// to the next block in the block chain.
		if ( lastelement == lastblockend )
		{
			// If there are not free blocks then add a new block to the chain.
			if ( lastblock->next == firstblock )
				lastblock->next = new Block( lastblock->next );

			lastblock = lastblock->next;
			lastblockend = (void*)( (UInt)lastblock + sizeof(Block) );
			lastelement = &lastblock->data[0];
		};
	};


	bool Empty() const
	{
		return firstelement == lastelement;
	};


	void Pop(DataType& data)
	{
		Assert( not Empty() );

		data = *firstelement;
		firstelement++;

		// If we hit the end of the block then we need to move the pointers
		// to the next block in the block chain.
		if ( firstelement == firstblockend )
		{
			firstblock = firstblock->next;
			firstblockend = (void*)( (UInt)firstblock + sizeof(Block) );
			firstelement = &firstblock->data[0];

			// If there is more than one free block then delete the excess block.
			if ( lastblock->next->next != firstblock )
			{
				Block* blocktodelete = lastblock->next;
				lastblock->next = lastblock->next->next;
				delete blocktodelete;
			};
		};
	};


	DataType Pop()
	{
		DataType data;
		Pop(data);
		return data;
	};


private:

	Queue(const Queue&) {};

	class Block
	{
	public:

		Block* next;
		DataType data[blocksize];

		Block()
		{
			next = NULL;
		};
		
		Block(Block* nextblock)
		{
			next = nextblock;
		};
	};


	Block* firstblock;
	void* firstblockend;
	Block* lastblock;
	void* lastblockend;
	DataType* firstelement;
	DataType* lastelement;
};


} // Buffers
} // dHLT

#endif // dHLT_BUFFERS_QUEUE_HPP
