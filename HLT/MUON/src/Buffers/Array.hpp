////////////////////////////////////////////////////////////////////////////////
//
// Author: Artur Szostak
// Email:  artur@alice.phy.uct.ac.za | artursz@iafrica.com
//
////////////////////////////////////////////////////////////////////////////////

#ifndef dHLT_BUFFERS_ARRAY_HPP
#define dHLT_BUFFERS_ARRAY_HPP

#include "BasicTypes.hpp"
#include "Utils.hpp"
#include "new.hpp"

namespace dHLT
{
namespace Buffers
{


template <typename DataType>
class Array
{
protected:

	UInt size;
	DataType* array;

public:

	class Iterator
	{
	public:

		Iterator() {};

		Iterator(const Iterator& iter)
		{
			it = iter.it;
		};

		Iterator(DataType* iter)
		{
			it = iter;
		};

		DataType& operator * ()
		{
			return *it;
		};

		const DataType& operator * () const
		{
			return *it;
		};

		DataType* operator -> ()
		{
			return it;
		};

		const DataType* operator -> () const
		{
			return it;
		};

		Iterator& operator ++ ()
		{
			it++;
			return *this;
		};

		Iterator operator ++ (int)
		{
			Iterator copy = *this;
			it++;
			return copy;
		};

		Iterator& operator -- ()
		{
			it--;
			return *this;
		};

		Iterator operator -- (int)
		{
			Iterator copy = *this;
			it--;
			return copy;
		};

		operator DataType* ()
		{
			return it;
		};

		operator const DataType* () const
		{
			return it;
		};

		friend bool operator == (const Iterator& a, const Iterator& b)
		{
			return a.it == b.it;
		};

		friend bool operator != (const Iterator& a, const Iterator& b)
		{
			return a.it != b.it;
		};


	protected:

		DataType* it;
	};

	
	Array()
	{
		size = 0;
		array = NULL;
	};
	
	
	~Array()
	{
		if (array != NULL)
			delete [] array;
	};
	

	bool Empty() const
	{
		return Size() == 0;
	};


	void Add(const DataType& data)
	{
		DataType* tmp = new DataType[size+1];
		for (UInt i = 0; i < size; i++)
			tmp[i] = array[i];
		tmp[size++] = data;
		if (array != NULL) delete [] array;
		array = tmp;
	};


	DataType* New()
	{
		DataType* tmp = new DataType[size+1];
		for (UInt i = 0; i < size; i++)
			tmp[i] = array[i];
		size++;
		if (array != NULL) delete [] array;
		array = tmp;
		return &array[size-1];
	};
	
	
	inline DataType* Add()
	{
		return New();
	};
	

	DataType* AddNew(const DataType& data)
	{
		DataType* tmp = new DataType[size+1];
		for (UInt i = 0; i < size; i++)
			tmp[i] = array[i];
		tmp[size++] = data;
		if (array != NULL) delete [] array;
		array = tmp;
		return &array[size-1];
	};


	DataType* AddUniquely(const DataType& data)
	{
		DataType* result = Find(data);
		if (result == NULL)
			return Add(data);
		else
			return result;
	};


	void Remove(const DataType& data)
	{
		UInt index = Size();
		for (UInt i = 0; i < Size(); i++)
		{
			if ( array[i] == data )
			{
				index = i;
				break;
			};
		};
		if (index < Size())
		{
			DataType* tmp = new DataType[size-1];
			UInt i;
			for (i = 0; i < index; i++)
				tmp[i] = array[i];
			for (i = index + 1; i < Size(); i++)
				tmp[i-1] = array[i];
			if (array != NULL) delete [] array;
			array = tmp;
			size--;
		};
	};


	DataType* Find(const DataType& data)
	{
		for (UInt i = 0; i < Size(); i++)
		{
			if ( array[i] == data )
				return &array[i];
		};
		return NULL;
	};


	const DataType* Find(const DataType& data) const
	{
		for (UInt i = 0; i < Size(); i++)
		{
			if ( array[i] == data )
				return &array[i];
		};
		return NULL;
	};


	bool Contains(const DataType& data) const
	{
		return Find(data) != NULL;
	};


	void Clear()
	{
		size = 0;
		if (array != NULL)
		{
			delete [] array;
			array = NULL;
		};
	};


	UInt Size() const
	{
		return size;
	};


	UInt Count() const
	{
		return size;
	};


	Iterator First()
	{
		return &array[0];
	};


	const Iterator First() const
	{
		return &array[0];
	};


	Iterator End()
	{
		return NULL;
	};


	const Iterator End() const
	{
		return NULL;
	};
	
	
	DataType& operator [] (const UInt index)
	{
		Assert( index < size );
		return array[index];
	};
	
	
	const DataType& operator [] (const UInt index) const
	{
		Assert( index < size );
		return array[index];
	};

};


} // Buffers
} // dHLT

#endif // dHLT_BUFFERS_ARRAY_HPP
