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
		return QueueAsyncTask(&QueueAsyncMemberTaskHelper, new AliHLTAsyncMemberProcessorContainer(obj, function, data));
	}

	void* InitializeAsyncMemberTask(T* obj, void* (T::*function)(void*), void* data)
	{
		return InitializeAsyncTask(&QueueAsyncMemberTaskHelper, new AliHLTAsyncMemberProcessorContainer(obj, function, data));
	}

private:
	AliHLTAsyncMemberProcessor(const AliHLTAsyncMemberProcessor&);
	AliHLTAsyncMemberProcessor& operator=(const AliHLTAsyncMemberProcessor&);

	struct AliHLTAsyncMemberProcessorContainer
	{
		AliHLTAsyncMemberProcessorContainer(T* obj, void* (T::*function)(void*), void* data) : fObj(obj), fFunction(function), fData(data) {}
		T* fObj;
		void* (T::*fFunction)(void*);
		void* fData;
	};

	static void* QueueAsyncMemberTaskHelper(void* data)
	{
		AliHLTAsyncMemberProcessorContainer *tmpData = (AliHLTAsyncMemberProcessorContainer*) data;
		AliHLTAsyncMemberProcessorContainer tmpDataCopy = *tmpData;
		delete tmpData;
		return((tmpDataCopy.fObj->*tmpDataCopy.fFunction)(tmpDataCopy.fData));
	}
};

#endif
