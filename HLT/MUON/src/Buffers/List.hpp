////////////////////////////////////////////////////////////////////////////////
//
// Author: Artur Szostak
// Email:  artur@alice.phy.uct.ac.za | artursz@iafrica.com
//
////////////////////////////////////////////////////////////////////////////////

#ifndef dHLT_BUFFERS_LIST_HPP
#define dHLT_BUFFERS_LIST_HPP

#include "BasicTypes.hpp"
#include "Utils.hpp"
#include "new.hpp"

#ifdef DEBUG
#include <iostream>
#include "Debug/print.hpp"
using std::endl;
using std::cout;
#endif // DEBUG

namespace dHLT
{
namespace Buffers
{


template <typename DataType>
class List
{
protected:

	struct Node
	{
		Node* next;
		DataType data;
	};

	Node* first;

public:

	class ConstIterator
	{
	public:

		ConstIterator()
		{
			current = NULL;
		};

		ConstIterator(const ConstIterator& iter)
		{
			current = iter.current;
		};

		ConstIterator(Node* node)
		{
			current = node;
		};

		const DataType& operator * () const
		{
			return current->data;
		};

		const DataType* operator -> () const
		{
			return &current->data;
		};

		ConstIterator& operator ++ ()
		{
			Assert( current != NULL );
			current = current->next;
			return *this;
		};

		ConstIterator operator ++ (int)
		{
			Assert( current != NULL );
			ConstIterator copy = *this;
			current = current->next;
			return copy;
		};
		
		operator const DataType* () const
		{
			return &current->data;
		};

		friend bool operator == (const ConstIterator& a, const ConstIterator& b)
		{
			return a.current == b.current;
		};

		friend bool operator != (const ConstIterator& a, const ConstIterator& b)
		{
			return a.current != b.current;
		};

	protected:
	
		friend class List;

		Node* current;
	};
	

	class Iterator : public ConstIterator
	{
	public:

		Iterator() : ConstIterator()
		{
			previous = NULL;
		};

		Iterator(const Iterator& iter) : ConstIterator(iter)
		{
			previous = iter.previous;
		};

		Iterator(Node* current, Node* prev) : ConstIterator(current)
		{
			previous = prev;
		};

		DataType& operator * ()
		{
			return ConstIterator::current->data;
		};

		DataType* operator -> ()
		{
			return &ConstIterator::current->data;
		};

		Iterator& operator ++ ()
		{
			Assert( ConstIterator::current != NULL );
			previous = ConstIterator::current;
			ConstIterator::current = ConstIterator::current->next;
			return *this;
		};

		Iterator operator ++ (int)
		{
			Assert( ConstIterator::current != NULL );
			Iterator copy = *this;
			previous = ConstIterator::current;
			ConstIterator::current = ConstIterator::current->next;
			return copy;
		};
		
		operator DataType* ()
		{
			if (ConstIterator::current != NULL)
				return &ConstIterator::current->data;
			else
				return NULL;
		};

	protected:
	
		friend class List;

		Node* previous;
	};

	
	List()
	{
		first = NULL;
	};
	
	
	~List()
	{
		Node* current = first;
		while (current != NULL)
		{
			Node* temp = current;
			current = current->next;
			delete temp;
		};
	};
	

	bool Empty() const
	{
		return first == NULL;
	};

	
	DataType* New()
	{
		Node* newnode = new Node();
		newnode->next = first;
		first = newnode;
		return &newnode->data;
	};
	
	
	inline DataType* Add()
	{
		return New();
	};
	

	void Add(const DataType& data)
	{
		DataType* newdata = New();
		*newdata = data;
	};
	

	DataType* AddNew(const DataType& data)
	{
		DataType* newdata = Add();
		*newdata = data;
		return newdata;
	};


	DataType* AddUniquely(const DataType& data)
	{
		DataType* result = Find(data);
		if (result == NULL)
			return AddNew(data);
		else
			return result;
	};

	
	void Remove(const UInt index)
	{
		Node* previous = NULL;
		Node* current = first;
		for (UInt i = 0; i < index; i++)
		{
			Assert( current != NULL );
			previous = current;
			current = current->next;
		};
		Node* temp;
		if (previous == NULL)
		{
			temp = first;
			first = first->next;
		}
		else
		{
			temp = current;
			previous->next = current->next;
		};
		delete temp;
	};
	

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
	};
	

	void Remove(Iterator& iter)
	{
		Node* previous = iter.previous;
		Node* current = iter.current;
		Assert( current != NULL );
		
		Node* temp;
		if (previous == NULL)
		{
			temp = first;
			first = first->next;
		}
		else
		{
			temp = current;
			previous->next = current->next;
		};
		delete temp;
	};


	Iterator Find(const DataType& data)
	{
		Node* previous = NULL;
		Node* current = first;
		while (current != NULL)
		{
			if (current->data == data)
				return Iterator(current, previous);
			previous = current;
			current = current->next;
		};
		return End();
	};


	ConstIterator Find(const DataType& data) const
	{
		Node* current = first;
		while (current != NULL)
		{
			if (current->data == data)
				return current;
			current = current->next;
		};
		return End();
	};


	template <typename PredicateType>
	Iterator Find(PredicateType predicate)
	{
		Node* previous = NULL;
		Node* current = first;
		while (current != NULL)
		{
			if ( predicate(current->data) )
				return Iterator(current, previous);
			previous = current;
			current = current->next;
		};
		return End();
	};


	template <typename PredicateType>
	ConstIterator Find(PredicateType predicate) const
	{
		Node* current = first;
		while (current != NULL)
		{
			if ( predicate(current->data) )
				return current;
			current = current->next;
		};
		return End();
	};


	bool Contains(const DataType& data) const
	{
		return Find(data) != End();
	};


	void Clear()
	{
		Node* current = first;
		while (current != NULL)
		{
			Node* temp = current;
			current = current->next;
			delete temp;
		};
		first = NULL;
	};


	Iterator First()
	{
		return Iterator(first, NULL);
	};


	ConstIterator First() const
	{
		return first;
	};


	Iterator End()
	{
		return Iterator(NULL, NULL);
	};


	ConstIterator End() const
	{
		return NULL;
	};

	
#	ifdef DEBUG

	void Dump()
	{
		Node* current = first;
		while (current != NULL)
		{
			cout << current->data << endl;
			current = current->next;
		};
	};

#	endif // DEBUG
};


} // Buffers
} // dHLT

#endif // dHLT_BUFFERS_LIST_HPP
