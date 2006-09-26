////////////////////////////////////////////////////////////////////////////////
//
// Author: Artur Szostak
// Email:  artur@alice.phy.uct.ac.za | artursz@iafrica.com
//
////////////////////////////////////////////////////////////////////////////////

#ifndef ALIHLTMUONCORECOUNTEDLIST_H
#define ALIHLTMUONCORECOUNTEDLIST_H

#include "AliHLTMUONBasicTypes.h"
#include "AliHLTMUONUtils.h"
#include "AliHLTMUONOutOfMemory.h"
#include "Buffers/List.hpp"

#ifdef DEBUG
#include <iostream>
#include "Debug/print.hpp"
using std::endl;
using std::cout;
#endif // DEBUG


template <typename DataType>
class AliHLTMUONCoreCountedList : public AliHLTMUONCoreList<DataType>
{
public:

	typedef typename AliHLTMUONCoreList<DataType>::Iterator Iterator;
	typedef typename AliHLTMUONCoreList<DataType>::ConstIterator ConstIterator;


	AliHLTMUONCoreCountedList() : AliHLTMUONCoreList<DataType>(), fCount(0)
	{
		fCount = 0;
	};

	
	DataType* New()
	{
		DataType* newdata = AliHLTMUONCoreList<DataType>::New();
		fCount++;
		return newdata;
	};
	
	
	DataType* Add()
	{
		return New();
	};
	

	void Add(const DataType& data)
	{
		AliHLTMUONCoreList<DataType>::Add(data);
		fCount++;
	};
	

	DataType* AddNew(const DataType& data)
	{
		DataType* newdata = AliHLTMUONCoreList<DataType>::AddNew(data);
		fCount++;
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
		AliHLTMUONCoreList<DataType>::Remove(index);
		fCount--;
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
		AliHLTMUONCoreList<DataType>::Remove(iter);
		fCount--;
	};


	void Clear()
	{
		AliHLTMUONCoreList<DataType>::Clear();
		fCount = 0;
	};
	
	
	UInt Count() const
	{
		return fCount;
	};


protected:

	UInt fCount;

private:

	// Hide copy constructor and assignment operator
	
	AliHLTMUONCoreCountedList(const AliHLTMUONCoreCountedList& list)
		: AliHLTMUONCoreList<DataType>(), fCount(0)
	{}

	
	AliHLTMUONCoreCountedList& operator = (const AliHLTMUONCoreCountedList& list)
	{
		return *this;
	}
};


#endif // ALIHLTMUONCORECOUNTEDLIST_H
