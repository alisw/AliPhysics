#ifndef ALIHLTMUONLIST_H
#define ALIHLTMUONLIST_H
/**************************************************************************
 * Copyright(c) 1998-2007, ALICE Experiment at CERN, All rights reserved. *
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

// $Id$

/**
 * @file   AliHLTMUONList.h
 * @author Artur Szostak <artursz@iafrica.com>
 * @date   29 May 2007
 * @brief  Declaration of a singly linked-list class which preallocates memory to give faster online performance.
 */

#include "AliHLTDataTypes.h"
#include <ostream>
#include <cassert>
#include <cstdlib>
#include <new>

/**
 * The AliHLTMUONList class is a unidirectional linked-list which preallocates
 * a large memory buffer where it stores its elements. This strategy gives higher
 * online performance.
 */
template <typename DataType>
class AliHLTMUONList
{
protected:

	struct Node;  // forward declaration for ConstIterator and Iterator

public:

	/**
	 * This is an iterator class which allows one to iterate through the list
	 * using the prefix or postfix operators ++.
	 */
	class ConstIterator
	{
	public:

		ConstIterator() : fCurrent(NULL) {}

		ConstIterator(const ConstIterator& iter) : fCurrent(iter.fCurrent) {}

		ConstIterator(Node* node) : fCurrent(node) {}

		ConstIterator& operator = (const ConstIterator& iter)
		{
		        if (this==&iter) return *this;
			fCurrent = iter.fCurrent;
			return *this;
		}
		
		virtual ~ConstIterator() {} // Just to make gcc -Weffc++ option shutup.

		// Dereference operator returns the elements data.
		const DataType& operator * () const { return fCurrent->fData; }
		
		// Pointer dereferencing returns the elements data.
		const DataType* operator -> () const { return &fCurrent->fData; }

		// Move to the next element in the list.
		ConstIterator& operator ++ ()
		{
			assert( fCurrent != NULL );
			fCurrent = fCurrent->fNext;
			return *this;
		}

		// Move to the next element in the list.
		ConstIterator operator ++ (int)
		{
			assert( fCurrent != NULL );
			ConstIterator copy = *this;
			fCurrent = fCurrent->fNext;
			return copy;
		}
		
		// Typecast operator returns a pointer to the data.
		operator const DataType* () const { return &fCurrent->fData; }

		friend bool operator == (const ConstIterator& a, const ConstIterator& b)
		{
			return a.fCurrent == b.fCurrent;
		}

		friend bool operator != (const ConstIterator& a, const ConstIterator& b)
		{
			return a.fCurrent != b.fCurrent;
		}

	protected:
	
		friend class AliHLTMUONList;
		Node* fCurrent; // The pointer to the current element selected by the iterator.
	};
	
	/**
	 * This is an iterator class which allows one to iterate through the list
	 * just like ConstIterator, but also allows modification of the elements.
	 */
	class Iterator : public ConstIterator
	{
	public:

		Iterator() : ConstIterator(), fPrevious(NULL) {}

		Iterator(const Iterator& iter) : ConstIterator(iter), fPrevious(iter.fPrevious) {}

		Iterator(Node* current, Node* prev) : ConstIterator(current), fPrevious(prev) {}

		Iterator& operator = (const Iterator& iter)
		{
		        if (this==&iter) return *this;
			ConstIterator::operator = (iter);
			fPrevious = iter.fPrevious;
			return *this;
		}
		
		virtual ~Iterator() {} // Just to make gcc -Weffc++ option shutup.

		// Dereference operator returns the elements data.
		DataType& operator * () { return ConstIterator::fCurrent->fData; }

		// Pointer dereferencing returns the elements data.
		DataType* operator -> () { return &ConstIterator::fCurrent->fData; }

		// Move to the next element in the list.
		Iterator& operator ++ ()
		{
			assert( ConstIterator::fCurrent != NULL );
			fPrevious = ConstIterator::fCurrent;
			ConstIterator::fCurrent = ConstIterator::fCurrent->fNext;
			return *this;
		}

		// Move to the next element in the list.
		Iterator operator ++ (int)
		{
			assert( ConstIterator::fCurrent != NULL );
			Iterator copy = *this;
			fPrevious = ConstIterator::fCurrent;
			ConstIterator::fCurrent = ConstIterator::fCurrent->fNext;
			return copy;
		}
		
		// Typecast operator returns a pointer to the data.
		operator DataType* ()
		{
			if (ConstIterator::fCurrent != NULL)
				return &ConstIterator::fCurrent->fData;
			else
				return NULL;
		}

	protected:
	
		friend class AliHLTMUONList;
		Node* fPrevious;  // The pointer to the previous element.
	};

	/**
	 * Constructor allocates a buffer to hold at least 'minentries' number
	 * of nodes for the list.
	 */
	AliHLTMUONList(AliHLTUInt32_t maxentries = 1024*4) :
		fFirst(NULL), fNextFree(0), fMaxEntries(maxentries), fEntry(NULL)
	{
		// Allocate buffer space.
		fEntry = reinterpret_cast<NodeEntry*>(
				malloc( sizeof(NodeEntry) * fMaxEntries )
			);
		if (fEntry == NULL) throw std::bad_alloc();
		// Set the availability flags.
		for (AliHLTUInt32_t i = 0; i < fMaxEntries; i++)
			fEntry[i].fFree = true;
	}
	
	/**
	 * Deep copy the list.
	 */
	AliHLTMUONList(const AliHLTMUONList& list) :
		fFirst(NULL), fNextFree(0), fMaxEntries(list.fMaxEntries),
		fEntry(NULL)
	{
		if (list.fFirst != NULL)
		{
			// First allocate the buffer space.
			fEntry = reinterpret_cast<NodeEntry*>(
					malloc( sizeof(NodeEntry) * fMaxEntries )
				);
			if (fEntry == NULL) throw std::bad_alloc();
			// Set the availability flags.
			for (AliHLTUInt32_t i = 0; i < fMaxEntries; i++)
				fEntry[i].fFree = true;
		
			// Now copy the list by adding new node entries.
			Add(list.fFirst->fData);
			Node* current = fFirst;
			Node* listcurrent = list.fFirst->fNext;
			while (listcurrent != NULL)
			{
				current->fNext = NewNode(listcurrent->fData);
				current = current->fNext;
				listcurrent = listcurrent->fNext;
			}
		}
	}
	
	/**
	 * Delete the list and release all allocated memory.
	 */
	virtual ~AliHLTMUONList()
	{
		Node* current = fFirst;
		while (current != NULL)
		{
			Node* temp = current;
			current = current->fNext;
			DeleteNode(temp);
		}
		// Delete the buffer space after all nodes are deleted,
		// else we will cause a crash.
		assert(fEntry != NULL);
		free(fEntry);
	}
	
	/**
	 * Deep copy the list.
	 */
	AliHLTMUONList& operator = (const AliHLTMUONList& list)
	{
		if (list.fFirst == NULL)
		{
			Clear();
			return *this;
		}
		if (fFirst == NULL)
		{
			// The list in 'this' object is empty so we need to create
			// new nodes for everything.
			Add(list.fFirst->fData);
			Node* current = fFirst;
			Node* listcurrent = list.fFirst->fNext;
			while (listcurrent != NULL)
			{
				current->fNext = NewNode(listcurrent->fData);
				current = current->fNext;
				listcurrent = listcurrent->fNext;
			}
			return *this;
		}
		
		// At least the first node does not need to be created.
		fFirst->fData = list.fFirst->fData;
		
		Node* current = fFirst;
		Node* listcurrent = list.fFirst->fNext;
		while (listcurrent != NULL)
		{
			if (current->fNext == NULL)
			{
				// The list of 'this' object is shorter so we need
				// to create new nodes.
				do
				{
					current->fNext = NewNode(listcurrent->fData);
					current = current->fNext;
					listcurrent = listcurrent->fNext;
				}
				while (listcurrent != NULL);
				break;
			}
			current->fNext->fData = listcurrent->fData;
			current = current->fNext;
			listcurrent = listcurrent->fNext;
		}
		
		// Unlink the remaining nodes.
		Node* temp = current->fNext;
		current->fNext = NULL;
		current = temp;
		// Remove any excess nodes from 'this' list.
		while (current != NULL)
		{
			temp = current;
			current = current->fNext;
			DeleteNode(temp);
		}
		
		return *this;
	}
	
	/**
	 * Returns true if the list is empty and false otherwise.
	 */
	bool Empty() const { return fFirst == NULL; }
	
	/**
	 * Adds a new element to the start of the linked list.
	 * @return  The pointer to the new element to fill its fields.
	 */
	DataType* Add()
	{
		Node* newnode = NewNode();
		newnode->fNext = fFirst;
		fFirst = newnode;
		return &newnode->fData;
	}
	
	/**
	 * Adds a new element to the start of the linked list and fills it with
	 * the data specified in 'data'.
	 * @param data  The value to set the new element to.
	 */
	void Add(const DataType& data)
	{
		Node* newnode = NewNode(data);
		newnode->fNext = fFirst;
		fFirst = newnode;
	}

	/**
	 * Searches the list if the element 'data' is already in the list. If it
	 * is then a pointer to the existing element is returned, otherwise a new
	 * element is created and a pointer to it is returned.
	 * @param data  The value to search for or set the new element to.
	 * @return  A pointer to the existing or new element.
	 */
	DataType* AddUniquely(const DataType& data)
	{
		Iterator result = Find(data);
		if (result == ConstIterator(NULL))
		{
			DataType* newdata = Add();
			*newdata = data;
			return newdata;
		}
		else
			return result;
	}

	/**
	 * Removes the index'th element from the list.
	 * No error checking is done so there better be at least 'index' number
	 * of entries in the list. You can use Count() to find out how many
	 * entries there are.
	 */
	void Remove(const AliHLTUInt32_t index)
	{
		Node* previous = NULL;
		Node* current = fFirst;
		for (AliHLTUInt32_t i = 0; i < index; i++)
		{
			assert( current != NULL );
			previous = current;
			current = current->fNext;
		}
		Node* temp;
		if (previous == NULL)
		{
			temp = fFirst;
			fFirst = fFirst->fNext;
		}
		else
		{
			temp = current;
			previous->fNext = current->fNext;
		}
		DeleteNode(temp);
	}
	
	/**
	 * Looks for the entry with the same values as 'data' and removes it
	 * from the list. If the entry could not be found then false is returned.
	 * However if it is found then it is deleted and true is returned.
	 */
	bool Remove(const DataType& data)
	{
		Iterator current = Find(data);
		if (current != ConstIterator(NULL))
		{
			Remove(current);
			return true;
		}
		else
			return false;
	}
	
	/**
	 * Removes the entry pointed to by the iterator which must have been
	 * extracted from the list with a call to First() and/or several calls
	 * to the iterators increment operators.
	 * @param iter  The entry to remove from the list.
	 */
	void Remove(Iterator& iter)
	{
		Node* previous = iter.fPrevious;
		Node* current = iter.fCurrent;
		assert( current != NULL );
		
		Node* temp;
		if (previous == NULL)
		{
			temp = fFirst;
			fFirst = fFirst->fNext;
		}
		else
		{
			temp = current;
			previous->fNext = current->fNext;
		}
		DeleteNode(temp);
	}

	/**
	 * Finds the first element with the same value as 'data' and returns an
	 * iterator to it. The iterator points to the end of the list i.e. End()
	 * if nothing was found.
	 * @param data  The value to look for.
	 * @return  The iterator pointing to the found element or End().
	 */
	Iterator Find(const DataType& data)
	{
		Node* previous = NULL;
		Node* current = fFirst;
		while (current != NULL)
		{
			if (current->fData == data)
				return Iterator(current, previous);
			previous = current;
			current = current->fNext;
		}
		return End();
	}

	ConstIterator Find(const DataType& data) const
	{
		Node* current = fFirst;
		while (current != NULL)
		{
			if (current->fData == data)
				return ConstIterator(current);
			current = current->fNext;
		};
		return End();
	}

	/**
	 * Finds the first element in the list that returns true for the predicate.
	 * That is, the first element for which the data evaluates to true in the
	 * test: predicate(fData). The iterator points to the end of the list
	 * i.e. End() if nothing was found.
	 * @param predicate  Predicate that the data must return true for.
	 *                   This can usually be some one variable function.
	 * @return  The iterator pointing to the found element or End().
	 */
	template <typename PredicateType>
	Iterator Find(PredicateType predicate)
	{
		Node* previous = NULL;
		Node* current = fFirst;
		while (current != NULL)
		{
			if ( predicate(current->fData) )
				return Iterator(current, previous);
			previous = current;
			current = current->fNext;
		}
		return End();
	}

	template <typename PredicateType>
	ConstIterator Find(PredicateType predicate) const
	{
		Node* current = fFirst;
		while (current != NULL)
		{
			if ( predicate(current->fData) )
				return ConstIterator(current);
			current = current->fNext;
		}
		return End();
	}

	/**
	 * Returns true if the list contains the specified element, else false.
	 * @param data  The value to search for in the list.
	 */
	bool Contains(const DataType& data) const { return Find(data) != End(); }

	/**
	 * This deletes all elements from the list.
	 */
	void Clear()
	{
		Node* current = fFirst;
		while (current != NULL)
		{
			Node* temp = current;
			current = current->fNext;
			DeleteNode(temp);
		}
		fFirst = NULL;
	}
	
	/**
	 * This deletes all elements from the list and resizes the buffer which
	 * is used to store the entries for the list.
	 */
	void Clear(AliHLTUInt32_t maxentries) throw(std::bad_alloc)
	{
		Clear();
		
		NodeEntry* newp = reinterpret_cast<NodeEntry*>(
				realloc(fEntry, sizeof(NodeEntry) * fMaxEntries)
			);
		if (newp == NULL) throw std::bad_alloc();
		
		// Set the availability flags for the node entries
		for (AliHLTUInt32_t i = fMaxEntries; i < maxentries; i++)
			fEntry[i].fFree = true;
	
		fEntry = newp;
		fMaxEntries = maxentries;
		fNextFree = 0;
	}

	/**
	 * Returns an iterator to the first element of the list.
	 */
	Iterator First() { return Iterator(fFirst, NULL); }
	ConstIterator First() const { return fFirst; }

	/**
	 * Returns an iterator pointing to the end of the list. The returned
	 * iterator does not actually point to a real data value but rather a
	 * sentinel value, and so it should not be dereferenced directly.
	 */
	Iterator End() { return Iterator(NULL, NULL); }
	ConstIterator End() const { return ConstIterator(NULL); }

	/**
	 * Counts and returns the number of elements in the list.
	 */
	AliHLTUInt32_t Count() const
	{
		AliHLTUInt32_t n = 0;
		Node* current = fFirst;
		while (current != NULL)
		{
			n++;
			current = current->fNext;
		}
		return n;
	}

protected:

	/**
	 * The internal node structure for the linked list.
	 */
	struct Node
	{
		Node* fNext;     // Next element in the list.
		DataType fData;  // The data for the current link.

		Node() : fNext(NULL), fData() {};
		Node(const DataType& data) : fNext(NULL), fData(data) {};
		Node(const Node& node) : fNext(node.fNext), fData(node.data) {};
		
		// Shallow copy the node.
		Node& operator = (const Node& node)
		{
			fNext = node.fNext;
			fData = node.fData;
			return *this;
		}
	};

	Node* fFirst;  // The first element in the linked list.
	
	struct NodeEntry
	{
		bool fFree; // Indicates if this block is free.
		Node fNode; // The node structure.
	};
	
	AliHLTUInt32_t fNextFree;   // The next entry that is presumably free.
	AliHLTUInt32_t fMaxEntries; // The number of node entries that can be stored in fEntries.
	NodeEntry* fEntry;          // Buffer of preallocated node entries.
	
	/**
	 * Locates the next free location in the fEntry buffer, creates a new node
	 * at that location and returns the pointer.
	 */
	Node* NewNode() throw(std::bad_alloc)
	{
		//return new Node();
		assert( fNextFree < fMaxEntries );
		AliHLTUInt32_t i = fNextFree;
		while (not fEntry[i].fFree and i < fMaxEntries) i++;
		if (fEntry[i].fFree)
		{
			fEntry[i].fFree = false;
			fNextFree = (i+1) % fMaxEntries;
			return new (&fEntry[i].fNode) Node();
		}
		i = 0;
		while (not fEntry[i].fFree and i < fNextFree) i++;
		if (fEntry[i].fFree)
		{
			fEntry[i].fFree = false;
			fNextFree = (i+1) % fMaxEntries;
			return new (&fEntry[i].fNode) Node();
		}
		throw std::bad_alloc();
	}
	
	Node* NewNode(const DataType& data) throw(std::bad_alloc)
	{
		//return new Node(data);
		assert( fNextFree < fMaxEntries );
		AliHLTUInt32_t i = fNextFree;
		while (not fEntry[i].fFree and i < fMaxEntries) i++;
		if (fEntry[i].fFree)
		{
			fEntry[i].fFree = false;
			fNextFree = (i+1) % fMaxEntries;
			return new (&fEntry[i].fNode) Node(data);
		}
		i = 0;
		while (not fEntry[i].fFree and i < fNextFree) i++;
		if (fEntry[i].fFree)
		{
			fEntry[i].fFree = false;
			fNextFree = (i+1) % fMaxEntries;
			return new (&fEntry[i].fNode) Node(data);
		}
		throw std::bad_alloc();
	}
	
	/**
	 * Destructs the node object.
	 */
	void DeleteNode(Node* node)
	{
		//delete node;
		node->~Node();
		// Now mark the entry as free.
		char* p = reinterpret_cast<char*>(node);
		p -= (sizeof(NodeEntry) - sizeof(Node));
		NodeEntry* entry = reinterpret_cast<NodeEntry*>(p);
		entry->fFree = true;
	}
};


/**
 * Stream operator to write the list to std::ostream. The output is printed in
 * the following format: [element1, element2, ..., elementN]
 */
template <typename DataType>
std::ostream& operator << (
		std::ostream& stream, const AliHLTMUONList<DataType>& list
	)
{
	typename AliHLTMUONList<DataType>::ConstIterator current;
	current = list.First();
	stream << "[";
	if (current == list.End())
	{
		stream << "]";
		return stream;
	}
	stream << (*current++);
	while (current != list.End())
		stream << ", " << (*current++);
	stream << "]";
	return stream;
}

#endif // ALIHLTMUONLIST_H
