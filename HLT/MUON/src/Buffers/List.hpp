////////////////////////////////////////////////////////////////////////////////
//
// Author: Artur Szostak
// Email:  artur@alice.phy.uct.ac.za | artursz@iafrica.com
//
////////////////////////////////////////////////////////////////////////////////

#ifndef ALIHLTMUONCOREBUFFERSLIST_H
#define ALIHLTMUONCOREBUFFERSLIST_H

#include "AliHLTMUONBasicTypes.h"
#include "AliHLTMUONUtils.h"
#include "AliHLTMUONOutOfMemory.h"

#ifdef DEBUG
#include <iostream>
#include "Debug/print.hpp"
using std::endl;
using std::cout;
#endif // DEBUG


template <typename DataType>
class AliHLTMUONCoreList
{
protected:

	struct Node
	{
		Node* fNext;
		DataType fData;

		Node() : fNext(NULL), fData() {};

		Node(const Node& node) : fNext(node.fNext), fData(node.data) {};

		Node& operator = (const Node& node)
		{
			fNext = node.fNext;
			fData = node.data;
			return *this;
		};
	};

	Node* fFirst;

public:


	class ConstIterator
	{
	public:

		ConstIterator() : fCurrent(NULL) {}

		ConstIterator(const ConstIterator& iter) : fCurrent(iter.fCurrent) {}

		ConstIterator(Node* node) : fCurrent(node) {}

		ConstIterator& operator = (const ConstIterator& iter)
		{
			fCurrent = iter.fCurrent;
			return *this;
		}

		const DataType& operator * () const
		{
			return fCurrent->fData;
		}

		const DataType* operator -> () const
		{
			return &fCurrent->fData;
		}

		ConstIterator& operator ++ ()
		{
			Assert( fCurrent != NULL );
			fCurrent = fCurrent->fNext;
			return *this;
		}

		ConstIterator operator ++ (int)
		{
			Assert( fCurrent != NULL );
			ConstIterator copy = *this;
			fCurrent = fCurrent->fNext;
			return copy;
		}
		
		operator const DataType* () const
		{
			return &fCurrent->fData;
		}

		friend bool operator == (const ConstIterator& a, const ConstIterator& b)
		{
			return a.fCurrent == b.fCurrent;
		}

		friend bool operator != (const ConstIterator& a, const ConstIterator& b)
		{
			return a.fCurrent != b.fCurrent;
		}

	protected:
	
		friend class AliHLTMUONCore;

		Node* fCurrent;
	};
	

	class Iterator : public ConstIterator
	{
	public:

		Iterator() : ConstIterator(), fPrevious(NULL) {}

		Iterator(const Iterator& iter) : ConstIterator(iter), fPrevious(iter.fPrevious) {}

		Iterator(Node* current, Node* prev) : ConstIterator(current), fPrevious(prev) {}

		Iterator& operator = (const Iterator& iter)
		{
			ConstIterator::operator = (iter);
			fPrevious = iter.fPrevious;
			return *this;
		}

		DataType& operator * ()
		{
			return ConstIterator::fCurrent->fData;
		}

		DataType* operator -> ()
		{
			return &ConstIterator::fCurrent->fData;
		}

		Iterator& operator ++ ()
		{
			Assert( ConstIterator::fCurrent != NULL );
			fPrevious = ConstIterator::fCurrent;
			ConstIterator::fCurrent = ConstIterator::fCurrent->fNext;
			return *this;
		}

		Iterator operator ++ (int)
		{
			Assert( ConstIterator::fCurrent != NULL );
			Iterator copy = *this;
			fPrevious = ConstIterator::fCurrent;
			ConstIterator::fCurrent = ConstIterator::fCurrent->fNext;
			return copy;
		}
		
		operator DataType* ()
		{
			if (ConstIterator::fCurrent != NULL)
				return &ConstIterator::fCurrent->fData;
			else
				return NULL;
		}

	protected:
	
		friend class AliHLTMUONCore;

		Node* fPrevious;
	};

	
	AliHLTMUONCoreList() : fFirst(NULL) {}
	
	
	~AliHLTMUONCoreList()
	{
		Node* current = fFirst;
		while (current != NULL)
		{
			Node* temp = current;
			current = current->fNext;
			delete temp;
		}
	}
	

	bool Empty() const
	{
		return fFirst == NULL;
	}

	
	DataType* New()
	{
		Node* newnode = new Node();
		newnode->fNext = fFirst;
		fFirst = newnode;
		return &newnode->fData;
	}
	
	
	DataType* Add()
	{
		return New();
	}
	

	void Add(const DataType& data)
	{
		DataType* newdata = New();
		*newdata = data;
	}
	

	DataType* AddNew(const DataType& data)
	{
		DataType* newdata = Add();
		*newdata = data;
		return newdata;
	}


	DataType* AddUniquely(const DataType& data)
	{
		DataType* result = Find(data);
		if (result == NULL)
			return AddNew(data);
		else
			return result;
	}

	
	void Remove(const UInt index)
	{
		Node* previous = NULL;
		Node* current = fFirst;
		for (UInt i = 0; i < index; i++)
		{
			Assert( current != NULL );
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
		delete temp;
	}
	

	bool Remove(const DataType& data)
	{
		Iterator current = Find(data);
		if ( current != ConstIterator(NULL) )
		{
			Remove(current);
			return true;
		}
		else
			return false;
	}
	

	void Remove(Iterator& iter)
	{
		Node* previous = iter.previous;
		Node* current = iter.current;
		Assert( current != NULL );
		
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
		delete temp;
	}


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
				return current;
			current = current->fNext;
		};
		return End();
	}


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
				return current;
			current = current->fNext;
		}
		return End();
	}


	bool Contains(const DataType& data) const
	{
		return Find(data) != End();
	}


	void Clear()
	{
		Node* current = fFirst;
		while (current != NULL)
		{
			Node* temp = current;
			current = current->fNext;
			delete temp;
		};
		fFirst = NULL;
	}


	Iterator First()
	{
		return Iterator(fFirst, NULL);
	}


	ConstIterator First() const
	{
		return fFirst;
	}


	Iterator End()
	{
		return Iterator(NULL, NULL);
	}


	ConstIterator End() const
	{
		return NULL;
	}

	
#	ifdef DEBUG

	void Dump()
	{
		Node* current = fFirst;
		while (current != NULL)
		{
			cout << current->fData << endl;
			current = current->fNext;
		}
	}

#	endif // DEBUG

private:

	// Hide copy contructor and assignment operator.
	AliHLTMUONCoreList(const AliHLTMUONCoreList& list) : fFirst(NULL) {}
	
	AliHLTMUONCoreList& operator = (const AliHLTMUONCoreList& list)
	{
		return *this;
	}

};


#endif // ALIHLTMUONCOREBUFFERSLIST_H
