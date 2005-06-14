////////////////////////////////////////////////////////////////////////////////
//
// Author: Artur Szostak
// Email:  artur@alice.phy.uct.ac.za | artursz@iafrica.com
//
////////////////////////////////////////////////////////////////////////////////

#ifndef dHLT_BUFFERS_COUNTED_LIST_HPP
#define dHLT_BUFFERS_COUNTED_LIST_HPP

#include "BasicTypes.hpp"
#include "Utils.hpp"
#include "new.hpp"
#include "Buffers/List.hpp"

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
class CountedList : public List<DataType>
{
public:

	typedef typename List<DataType>::Iterator Iterator;
	typedef typename List<DataType>::ConstIterator ConstIterator;


	CountedList() : List<DataType>()
	{
		count = 0;
	};

	
	DataType* New()
	{
		DataType* newdata = List<DataType>::New();
		count++;
		return newdata;
	};
	
	
	inline DataType* Add()
	{
		return New();
	};
	

	void Add(const DataType& data)
	{
		List<DataType>::Add(data);
		count++;
	};
	

	DataType* AddNew(const DataType& data)
	{
		DataType* newdata = List<DataType>::AddNew(data);
		count++;
		return newdata;
	};


	DataType* AddUniquely(const DataType& data)
	{
		DataType* result = Find(data);
		if (result == NULL)
			return Add(data);
		else
			return result;
	};

	
	void Remove(const UInt index)
	{
		List<DataType>::Remove(index);
		count--;
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
		List<DataType>::Remove(iter);
		count--;
	};


	void Clear()
	{
		List<DataType>::Clear();
		count = 0;
	};
	
	
	UInt Count() const
	{
		return count;
	};


protected:

	UInt count;
};


} // Buffers
} // dHLT

#endif // dHLT_BUFFERS_COUNTED_LIST_HPP
