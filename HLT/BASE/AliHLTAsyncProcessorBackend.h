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

//0: User Mutex, 1: Operation Mutex, 2: Input Mutex, 3: Output Mutex, 4: Output Full Mutex, 5: Buffer Mutex
#define ASYNC_MUTEX_COUNT 6

#define HLT_ASYNC_USE_SEM_T //Use POSIX semaphores instead of PTHREAD mutexes

#include <pthread.h>
#include <sys/wait.h>
#include <sys/types.h>
#include <unistd.h>
#include <errno.h>
#include <signal.h>
#ifdef HLT_ASYNC_USE_SEM_T
#include <semaphore.h>
#endif

class AliHLTAsyncProcessorBackend
{
public:
	AliHLTAsyncProcessorBackend() : fInitialized(false), fAsyncThread(0), fAsyncPID(0) {};
	~AliHLTAsyncProcessorBackend()
	{
		if (!fInitialized) return;
		for (int i = 0;i < ASYNC_MUTEX_COUNT;i++)
		{
#ifdef HLT_ASYNC_USE_SEM_T
			sem_destroy(&fMutexes[i]);
#else
			pthread_mutex_destroy(&fMutexes[i]);
#endif
		}
	};

	int Initialize(bool shared)
	{
#ifdef HLT_ASYNC_USE_SEM_T
		memset(fMutexes, 0, sizeof(fMutexes[0]) * ASYNC_MUTEX_COUNT);
#else
		pthread_mutexattr_t attr;
		pthread_mutexattr_init(&attr);
		if (shared) if (pthread_mutexattr_setpshared(&attr, PTHREAD_PROCESS_SHARED)) return(1);
#endif
		for (int i = 0;i < ASYNC_MUTEX_COUNT;i++)
		{
#ifdef HLT_ASYNC_USE_SEM_T
			if (sem_init(&fMutexes[i], shared, 1)) return(1);
#else
			if (pthread_mutex_init(&fMutexes[i], &attr)) return(1);
#endif
		}
#ifndef HLT_ASYNC_USE_SEM_T
		pthread_mutexattr_destroy(&attr);
#endif
		fInitialized = true;
		return(0);
	}

	int LockMutex(int i)
	{
#ifdef HLT_ASYNC_USE_SEM_T
		int retVal = sem_wait(&fMutexes[i]);
#else
		int retVal = pthread_mutex_lock(&fMutexes[i]);
#endif
		if (retVal)
		{
			fprintf(stderr, "Error locking mutex %d error %d\n", i, errno);
			usleep(10000);
		}
		return retVal;
	}

	int UnlockMutex(int i)
	{
#ifdef HLT_ASYNC_USE_SEM_T
		int retVal = sem_post(&fMutexes[i]);
#else
		int retVal = pthread_mutex_unlock(&fMutexes[i]);
#endif
		if (retVal)
		{
			fprintf(stderr, "Error unlocking mutex %d error %d\n", i, errno);
			usleep(100000);
		}
		return(retVal);
	}

	int TryLockMutex(int i)
	{
#ifdef HLT_ASYNC_USE_SEM_T
		return(sem_trywait(&fMutexes[i]));
#else
		return(pthread_mutex_trylock(&fMutexes[i]));
#endif
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
	
	int KillChildProcess(int waitTime)
	{
		do
		{
			int status;
			pid_t result = waitpid(fAsyncPID, &status, WNOHANG);
			if (result == 0)
			{
				usleep(1000);
			}
			else if (result == -1)
			{
				return(-1);
			}
			else
			{
				return(0);
			}
		} while (waitTime-- > 0);
		if (kill(fAsyncPID, SIGKILL)) return(-1);
		return(1);
	}

private:
#ifdef HLT_ASYNC_USE_SEM_T
	sem_t fMutexes[ASYNC_MUTEX_COUNT];
#else
	pthread_mutex_t fMutexes[ASYNC_MUTEX_COUNT];
#endif
	pthread_t fAsyncThread;
	pid_t fAsyncPID;
	bool fInitialized;
};

#undef ASYNC_MUTEX_COUNT
#ifdef HLT_ASYNC_USE_SEM_T
#undef HLT_ASYNC_USE_SEM_T
#endif

#endif
