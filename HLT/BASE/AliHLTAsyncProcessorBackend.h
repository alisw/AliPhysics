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
#include <sys/wait.h>
#include <sys/types.h>

class AliHLTAsyncProcessorBackend
{
public:
	AliHLTAsyncProcessorBackend() : fInitialized(false) {};
	~AliHLTAsyncProcessorBackend()
	{
		if (!fInitialized) return;
		for (int i = 0;i < ASYNC_MUTEX_COUNT;i++) pthread_mutex_destroy(&fMutexes[i]);
	};

	int Initialize(bool shared)
	{
		pthread_mutexattr_t attr;
		pthread_mutexattr_init(&attr);
		if (shared) if (pthread_mutexattr_setpshared(&attr, PTHREAD_PROCESS_SHARED)) return(1);
		for (int i = 0;i < ASYNC_MUTEX_COUNT;i++)
		{
			if (pthread_mutex_init(&fMutexes[i], &attr)) return(1);
		}
		pthread_mutexattr_destroy(&attr);

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

	int StartProcess(void* (*function)(void*), void* data)
	{
		int pid = fork();
		if (pid < 0) return(1);
		if (pid == 0)
		{
			function(data);
			exit(0);
		}
		fAsyncPID = pid;
		return(0);
	}
	
	int StopProcess()
	{
		int returnStatus;
		waitpid(fAsyncPID, &returnStatus, 0);
		return(returnStatus);
	}

private:
	bool fInitialized;
	pthread_mutex_t fMutexes[ASYNC_MUTEX_COUNT];
	pthread_t fAsyncThread;
	pid_t fAsyncPID;
};

#undef ASYNC_MUTEX_COUNT

#endif
