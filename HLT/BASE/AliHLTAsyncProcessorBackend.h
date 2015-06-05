#ifndef H_ALIHLTASYNCPROCESSORBACKEND
#define H_ALIHLTASYNCPROCESSORBACKEND

/* This file is property of and copyright by the ALICE HLT Project        * 
* ALICE Experiment at CERN, All rights reserved.                         *
* See cxx source for full Copyright notice                               */

/** @file    AliHLTAsyncProcessorBackend.h
@author  David Rohr (drohr@cern.ch)
*/

//This file is a backend for AliHLTAsyncProcessor, that implements the
//mutex and thread handling via pthreads

//0: User Mutex, 1: Operation Mutex, 2: Input Mutex, 3: Output Mutex, 4: Output Full Mutex
#define ASYNC_MUTEX_COUNT 5

#include <pthread.h>

class AliHLTAsyncProcessorBackend
{
public:
	AliHLTAsyncProcessorBackend() : fInitialized(false) {};
	~AliHLTAsyncProcessorBackend()
	{
		if (!fInitialized) return;
		for (int i = 0;i < ASYNC_MUTEX_COUNT;i++) pthread_mutex_destroy(&fMutexes[i]);
	};

	int Initialize()
	{
		for (int i = 0;i < ASYNC_MUTEX_COUNT;i++)
		{
			if (pthread_mutex_init(&fMutexes[i], NULL)) return(1);
		}

		fInitialized = true;
		return(0);
	}

	int LockMutex(int i)
	{
		return(pthread_mutex_lock(&fMutexes[i]));
	}

	int UnlockMutex(int i)
	{
		return(pthread_mutex_unlock(&fMutexes[i]));
	}

	int TryLockMutex(int i)
	{
		return(pthread_mutex_trylock(&fMutexes[i]));
	}

	int StartThread(void* (*function)(void*), void* data)
	{
		return(pthread_create(&fAsyncThread, NULL, function, data));
	}

	void* StopThread()
	{
		void* retVal;
		if (pthread_join(fAsyncThread, &retVal)) retVal = (void*) -1;
		return(retVal);
	}

private:
	bool fInitialized;
	pthread_mutex_t fMutexes[ASYNC_MUTEX_COUNT];
	pthread_t fAsyncThread;
};

#undef ASYNC_MUTEX_COUNT

#endif
