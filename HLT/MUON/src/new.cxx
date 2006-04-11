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
			fAllocList = NULL;
			fArrayAllocList = NULL;
		};
	
#ifdef CHECK_MEMORY_ALLOC_ON_EXIT
		~CheckMemoryAlloc()
		{
			Assert( fAllocList == NULL );
			Assert( fArrayAllocList == NULL );
		};
#endif // CHECK_MEMORY_ALLOC_ON_EXIT
		
		void AddAlloc(void* ptr)
		{
			Node* newnode = (Node*) malloc( sizeof(Node) );
			newnode->memory = ptr;
			newnode->next = fAllocList;
			fAllocList = newnode;
		};
		
		bool RemoveAlloc(void* ptr)
		{
			Node* previous = NULL;
			Node* current = fAllocList;
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
				fAllocList = current->next;
			free(( void*)current );
			return true;
		};
	
		void AddArrayAlloc(void* ptr)
		{
			Node* newnode = (Node*) malloc( sizeof(Node) );
			newnode->memory = ptr;
			newnode->next = fArrayAllocList;
			fArrayAllocList = newnode;
		};
		
		bool RemoveArrayAlloc(void* ptr)
		{
			Node* previous = NULL;
			Node* current = fArrayAllocList;
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
				fArrayAllocList = current->next;
			free( (void*)current );
			return true;
		};
		
	private:
	
		struct Node
		{
			void* memory;
			Node* next;
		};
	
		Node* fAllocList;
		Node* fArrayAllocList;
	
	} checkmem;
	
#endif // DEBUG

}; // end of namespace


void* operator new (size_t size) throw (std::bad_alloc)
{
	void* memory = malloc(size);
	if (memory == NULL) AliHLTMUONCoreThrowOutOfMemory();
	DebugMsg(99, "new(" << size << ") allocated: " << memory);
	DebugCode( checkmem.AddAlloc(memory) );
	return memory;
};


void* operator new [] (size_t size) throw (std::bad_alloc)
{
	void* memory = malloc(size);
	if (memory == NULL) AliHLTMUONCoreThrowOutOfMemory();
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

