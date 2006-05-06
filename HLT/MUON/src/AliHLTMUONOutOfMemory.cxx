/**************************************************************************
 * Copyright(c) 1998-1999, ALICE Experiment at CERN, All rights reserved. *
 *                                                                        *
 * Author: The ALICE Off-line Project.                                    *
 * Contributors are mentioned in the code where appropriate.              *
 *                                                                        *
 * Permission to use, copy, modify and distribute this software and its   *
 * documentation strictly for non-commercial purposes is hereby granted   *
 * without fee, provided that the above copyright notice appears in all   *
 * copies and that both the copyright notice and this permission notice   *
 * appear in the supporting documentation. The authors make no claims     *
 * about the suitability of this software for any purpose. It is          *
 * provided "as is" without express or implied warranty.                  *
 **************************************************************************/

/* $Id$ */

////////////////////////////////////////////////////////////////////////////////
//
// Author: Artur Szostak
// Email:  artur@alice.phy.uct.ac.za | artursz@iafrica.com
//
// AliHLTMUONOutOfMemory is defined to be used when the system runs
// out of memory. Do not throw this object directly but rather use the routine
// ThrowAliHLTMUONOutOfMemory which throws a pree allocated static object.
//
////////////////////////////////////////////////////////////////////////////////
 
#include "AliHLTMUONOutOfMemory.h"
#include "AliHLTMUONUtils.h"
#include <stdlib.h>


const char* AliHLTMUONOutOfMemory::Message() const throw()
{
	return "Out of memory.";
}

Int AliHLTMUONOutOfMemory::ErrorCode() const throw()
{
	return kAliHLTMUONOutOfMemory;
}

void ThrowAliHLTMUONOutOfMemory() throw (AliHLTMUONOutOfMemory)
{
	static AliHLTMUONOutOfMemory outofmemobject;
	throw outofmemobject;
}


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


void* operator new (size_t size) throw (AliHLTMUONOutOfMemory)
{
	void* memory = malloc(size);
	if (memory == NULL) ThrowAliHLTMUONOutOfMemory();
	DebugMsg(99, "new(" << size << ") allocated: " << memory);
	DebugCode( checkmem.AddAlloc(memory) );
	return memory;
};


void* operator new [] (size_t size) throw (AliHLTMUONOutOfMemory)
{
	void* memory = malloc(size);
	if (memory == NULL) ThrowAliHLTMUONOutOfMemory();
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

