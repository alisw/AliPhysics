#ifndef ALIHLTMUONCOUNTEDLIST_H
#define ALIHLTMUONCOUNTEDLIST_H
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
 * @file   AliHLTMUONCountedList.h
 * @author Artur Szostak <artursz@iafrica.com>
 * @date   
 * @brief  Declaration of a linked-list class which counts the number of
 *         elements it contains.
 */

#include "AliHLTMUONList.h"

/**
 * The AliHLTMUONCountedList class behaves just like the AliHLTMUONList class
 * but uses an internal counter to count the number of elements in the list.
 * This means calls to Count() are much more efficient.
 */
template <typename DataType>
class AliHLTMUONCountedList : public AliHLTMUONList<DataType>
{
public:

	typedef typename AliHLTMUONList<DataType>::Iterator Iterator;
	typedef typename AliHLTMUONList<DataType>::ConstIterator ConstIterator;
	
	AliHLTMUONCountedList(AliHLTUInt32_t maxentries = 1024*4) :
		AliHLTMUONList<DataType>(maxentries), fCount(0)
	{}

	// Perform a deep copy.	
	AliHLTMUONCountedList(const AliHLTMUONCountedList& list)
		: AliHLTMUONList<DataType>(list), fCount(list.fCount)
	{}

	// Perform a deep copy.
	AliHLTMUONCountedList& operator = (const AliHLTMUONCountedList& list)
	{
		AliHLTMUONList<DataType>::operator = (list);
		fCount = list.fCount;
		return *this;
	}
	
	virtual ~AliHLTMUONCountedList() {} // Just to make gcc -Weffc++ option shutup.
	
	/**
	 * Adds a new element to the start of the linked list.
	 * @return  The pointer to the new element to fill its fields.
	 */
	DataType* Add()
	{
		DataType* newdata = AliHLTMUONList<DataType>::Add();
		fCount++;
		return newdata;
	}

	/**
	 * Adds a new element to the start of the linked list and fills it with
	 * the data specified in 'data'.
	 * @param data  The value to set the new element to.
	 */
	void Add(const DataType& data)
	{
		AliHLTMUONList<DataType>::Add(data);
		fCount++;
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
		AliHLTMUONList<DataType>::Remove(index);
		fCount--;
	}
	
	/**
	 * Looks for the entry with the same values as 'data' and removes it
	 * from the list. If the entry could not be found then false is returned.
	 * However if it is found then it is deleted and true is returned.
	 */
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
	
	/**
	 * Removes the entry pointed to by the iterator which must have been
	 * extracted from the list with a call to First() and/or several calls
	 * to the iterators increment operators.
	 * @param iter  The entry to remove from the list.
	 */
	void Remove(Iterator& iter)
	{
		AliHLTMUONList<DataType>::Remove(iter);
		fCount--;
	}

	/**
	 * This deletes all elements from the list.
	 */
	void Clear()
	{
		AliHLTMUONList<DataType>::Clear();
		fCount = 0;
	}
	
	/**
	 * This deletes all elements from the list and resizes the buffer which
	 * is used to store the entries for the list.
	 */
	void Clear(AliHLTUInt32_t maxentries) throw(std::bad_alloc)
	{
		AliHLTMUONList<DataType>::Clear(maxentries);
		fCount = 0;
	}
	
	/**
	 * Counts and returns the number of elements in the list.
	 */
	AliHLTUInt32_t Count() const { return fCount; }

protected:

	AliHLTUInt32_t fCount;
};


#endif // ALIHLTMUONCOUNTEDLIST_H
