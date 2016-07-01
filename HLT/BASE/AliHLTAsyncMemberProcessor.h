#ifndef H_ALIHLTASYNCMEMBERPROCESSOR
#define H_ALIHLTASYNCMEMBERPROCESSOR

/* This file is property of and copyright by the ALICE HLT Project        * 
* ALICE Experiment at CERN, All rights reserved.                         *
* See cxx source for full Copyright notice                               */

/** @file    AliHLTAsyncMemberProcessor.h
@author  David Rohr (drohr@cern.ch)
*/

//This file is a wrapper for AliHLTAsyncProcessor, which allows to use
//member functions as callbacks in the queue.
//This class provides the two additional functions QueueAsyncMemberTask and
//InitializeAsyncMemberTask on top of the interface provided by AliHLTAsyncProcessor.
//They work in the same fashion as QueueAsyncTask and InitializeAsyncTask, but
//instead of non-member or static member function, they take member functions
//as callback paramters. Hence, the class instance (usually 'this') must be passed first.

#include "AliHLTAsyncProcessor.h"

template <class T>
class AliHLTAsyncMemberProcessor : public AliHLTAsyncProcessor
{
public:
	AliHLTAsyncMemberProcessor() : AliHLTAsyncProcessor() {}
	virtual ~AliHLTAsyncMemberProcessor() {}

	int QueueAsyncMemberTask(T* obj, void* (T::*function)(void*), void* data)
	{
		AliHLTAsyncMemberProcessorContainer* tmp = GetContainer();
		if (tmp == NULL) return(1);
		tmp->Set(obj, function, data);
		int retVal = QueueAsyncTask(&QueueAsyncMemberTaskHelper, tmp);
		if (retVal) FreeContainer(tmp);
		return (retVal);
	}

	int InitializeAsyncMemberTask(T* obj, void* (T::*function)(void*), void* data, void** pRetVal = NULL)
	{
		AliHLTAsyncMemberProcessorContainer* tmp = GetContainer();
		if (tmp == NULL) return(1);
		tmp->Set(obj, function, data);
		int retVal = InitializeAsyncTask(&QueueAsyncMemberTaskHelper, tmp, pRetVal);
		if (retVal) FreeContainer(tmp);
		return (retVal);
	}

private:
	AliHLTAsyncMemberProcessor(const AliHLTAsyncMemberProcessor&);
	AliHLTAsyncMemberProcessor& operator=(const AliHLTAsyncMemberProcessor&);
	
	struct AliHLTAsyncMemberProcessorContainer
	{
		void Set(T* obj, void* (T::*function)(void*), void* data) {fObj = obj; fFunction = function; fData = data; fUsed = true;}
		T* fObj;
		void* (T::*fFunction)(void*);
		void* fData;
		bool fUsed;
		bool fStaticallyAllocated;
	};

	AliHLTAsyncMemberProcessorContainer* GetContainer()
	{
		AliHLTAsyncMemberProcessorContainer* container;
		if (fMe->fQueueDepth == 0)
		{
			container = new AliHLTAsyncMemberProcessorContainer;
			container->fStaticallyAllocated = false;
			return(new AliHLTAsyncMemberProcessorContainer);
		}
		container = (AliHLTAsyncMemberProcessorContainer*) fMe->fChildBufferSpace;
		LockMutex(5);
		for (int i = 0;i <= fMe->fQueueDepth;i++)
		{
			if (!container[i].fUsed)
			{
				container[i].fUsed = true;
				container[i].fStaticallyAllocated = true;
				UnlockMutex(5);
				return(&container[i]);
			}
		}
		UnlockMutex(5);
		return(NULL);
	}
	
	static void FreeContainer(AliHLTAsyncMemberProcessorContainer* ptr)
	{
		if (!ptr->fStaticallyAllocated) delete ptr;
		else ptr->fUsed = false;
	}

	static void* QueueAsyncMemberTaskHelper(void* data)
	{
		AliHLTAsyncMemberProcessorContainer *tmpData = (AliHLTAsyncMemberProcessorContainer*) data;
		AliHLTAsyncMemberProcessorContainer tmpDataCopy = *tmpData;
		FreeContainer(tmpData);
		return((tmpDataCopy.fObj->*tmpDataCopy.fFunction)(tmpDataCopy.fData));
	}
	
	virtual size_t ChildSharedProcessBufferSize() {return sizeof(AliHLTAsyncMemberProcessorContainer) * (fMe->fQueueDepth + 1);}
};

#endif
