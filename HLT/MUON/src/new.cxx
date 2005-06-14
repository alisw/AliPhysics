////////////////////////////////////////////////////////////////////////////////
//
// Author: Artur Szostak
// Email:  artur@alice.phy.uct.ac.za | artursz@iafrica.com
//
////////////////////////////////////////////////////////////////////////////////

#include "new.hpp"
#include "Utils.hpp"
#include <stdlib.h>

namespace
{

#ifdef DEBUG

	/* Perform a check that on program termination all memory was released properly.
	   The internal list structures are implemented as linked lists using malloc and
	   free to allocate and release memory.
	   Can not use new or delete operators because thats what we are checking, if new
	   and delete were called properly.
	 */
	class CheckMemoryAlloc
	{
	public:
	
		CheckMemoryAlloc()
		{
			alloc_list = NULL;
			array_alloc_list = NULL;
		};
	
#ifdef CHECK_MEMORY_ALLOC_ON_EXIT
		~CheckMemoryAlloc()
		{
			Assert( alloc_list == NULL );
			Assert( array_alloc_list == NULL );
		};
#endif // CHECK_MEMORY_ALLOC_ON_EXIT
		
		void AddAlloc(void* ptr)
		{
			Node* newnode = (Node*) malloc( sizeof(Node) );
			newnode->memory = ptr;
			newnode->next = alloc_list;
			alloc_list = newnode;
		};
		
		bool RemoveAlloc(void* ptr)
		{
			Node* previous = NULL;
			Node* current = alloc_list;
			while (current != NULL)
			{
				if (current->memory == ptr) break;
				previous = current;
				current = current->next;
			};
			if (current == NULL) return false;
			if (previous != NULL)
				previous->next = current->next;
			else
				alloc_list = current->next;
			free(( void*)current );
			return true;
		};
	
		void AddArrayAlloc(void* ptr)
		{
			Node* newnode = (Node*) malloc( sizeof(Node) );
			newnode->memory = ptr;
			newnode->next = array_alloc_list;
			array_alloc_list = newnode;
		};
		
		bool RemoveArrayAlloc(void* ptr)
		{
			Node* previous = NULL;
			Node* current = array_alloc_list;
			while (current != NULL)
			{
				if (current->memory == ptr) break;
				previous = current;
				current = current->next;
			};
			if (current == NULL) return false;
			if (previous != NULL)
				previous->next = current->next;
			else
				array_alloc_list = current->next;
			free( (void*)current );
			return true;
		};
		
	private:
	
		struct Node
		{
			void* memory;
			Node* next;
		};
	
		Node* alloc_list;
		Node* array_alloc_list;
	
	} checkmem;
	
#endif // DEBUG

}; // end of namespace


void* operator new (size_t size) throw (dHLT::OutOfMemory)
{
	void* memory = malloc(size);
	if (memory == NULL) dHLT::ThrowOutOfMemory();
	DebugMsg(99, "new(" << size << ") allocated: " << memory);
	DebugCode( checkmem.AddAlloc(memory) );
	return memory;
};


void* operator new [] (size_t size) throw (dHLT::OutOfMemory)
{
	void* memory = malloc(size);
	if (memory == NULL) dHLT::ThrowOutOfMemory();
	DebugMsg(99, "new [] (" << size << ") allocated: " << memory);
	DebugCode( checkmem.AddArrayAlloc(memory) );
	return memory;
};


void operator delete (void* memory) throw ()
{
	DebugMsg(99, "delete(" << memory << ")");
	DebugCode(Assert( checkmem.RemoveAlloc(memory) ));
	free(memory);
};


void operator delete [] (void* memory) throw ()
{
	DebugMsg(99, "delete [] (" << memory << ")");
	DebugCode(Assert( checkmem.RemoveArrayAlloc(memory) ));
	free(memory);
};

