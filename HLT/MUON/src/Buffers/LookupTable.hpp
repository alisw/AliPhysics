////////////////////////////////////////////////////////////////////////////////
//
// Author: Artur Szostak
// Email:  artur@alice.phy.uct.ac.za | artursz@iafrica.com
//
////////////////////////////////////////////////////////////////////////////////

#ifndef dHLT_BUFFERS_LOOKUP_TABLE_HPP
#define dHLT_BUFFERS_LOOKUP_TABLE_HPP

#include "BasicTypes.hpp"
#include "Utils.hpp"
#include "new.hpp"

namespace dHLT
{
namespace Buffers
{


template <typename KeyType, class DataRecordType>
class LookupTable
{
private:

	struct Node
	{
		Node* next;
		DataRecordType data;

		Node(const KeyType key, Node* next = NULL)
			: data(key)
		{
			this->next = next;
		};
	};


	Node* first;


public:
	
	class Iterator
	{
	public:

		Iterator(Node* node = NULL)
		{
			current = node;
		};

		DataRecordType& operator * ()
		{
			return current->data;
		};

		const DataRecordType& operator * () const
		{
			return current->data;
		};

		DataRecordType* operator -> ()
		{
			return &current->data;
		};

		const DataRecordType* operator -> () const
		{
			return &current->data;
		};

		Iterator& operator ++ ()
		{
			current = current->next;
			return *this;
		};

		Iterator operator ++ (int)
		{
			Iterator copy = *this;
			current = current->next;
			return copy;
		};

		friend bool operator == (const Iterator& a, const Iterator& b)
		{
			return a.current == b.current;
		};

		friend bool operator != (const Iterator& a, const Iterator& b)
		{
			return a.current != b.current;
		};


	protected:

		Node* current;
	};


	LookupTable()
	{
		first = NULL;
	};

	~LookupTable()
	{
		Clear();
	};


	DataRecordType* Add(const KeyType& key)
	{
		first = new Node(key, first);
		return &first->data;
	};


	DataRecordType* AddUniquely(const KeyType& key)
	{
		DataRecordType* data = Find(key);
		if (data == NULL)
			return Add(key);
		else
			return data;
	};


	void Remove(const KeyType& key)
	{
		Node* previous = NULL;
		Node* current = first;
		if (current == NULL) return;
		while (current->data.Key() != key)
		{
			previous = current;
			current = current->next;
			if (current == NULL) return;
		};
		if (previous != NULL)
			previous->next = current->next;
		else
			first = current->next;
		delete current;
	};


	DataRecordType* Find(const KeyType& key)
	{
		Node* current = first;
		while (current != NULL)
		{
			if (current->data.Key() == key)
				return &current->data;
			current = current->next;
		};
		return NULL;
	};


	const DataRecordType* Find(const KeyType& key) const
	{
		Node* current = first;
		while (current != NULL)
		{
			if (current->data.Key() == key)
				return &current->data;
			current = current->next;
		};
		return NULL;
	};


	void Clear()
	{
		Node* current = first;
		while (current != NULL)
		{
			Node* nodetodelete = current;
			current = current->next;
			delete nodetodelete;
		};
		first = NULL;
	};


	Iterator First()
	{
		return Iterator(first);
	};

	const Iterator First() const
	{
		return Iterator(first);
	};

	Iterator End()
	{
		return NULL;
	};

	const Iterator End() const
	{
		return NULL;
	};
};


} // Buffers
} // dHLT

#endif // dHLT_BUFFERS_LOOKUP_TABLE_HPP
