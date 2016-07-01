#ifndef ALIHLTASYNCPROCESSOR_H
#define ALIHLTASYNCPROCESSOR_H

#define ALIHLTASYNCPROCESSOR_ALIGN 64

/* This file is property of and copyright by the ALICE HLT Project        * 
* ALICE Experiment at CERN, All rights reserved.                         *
* See cxx source for full Copyright notice                               */

/** @file    AliHLTAsyncProcessor.h
@author  David Rohr (drohr@cern.ch)
*/

//The AliHLTAsyncProcessor class implements an interface for asynchronous
//tasks for HLT components.

#include "AliHLTLogging.h"
class AliHLTAsyncProcessorBackend;
class AliHLTComponent;

class AliHLTAsyncProcessor : public AliHLTLogging
{
public:
	AliHLTAsyncProcessor();
	virtual ~AliHLTAsyncProcessor();

	//Initialize the Async Processor with given queue depth.
	//If initialize is not called, or depth is 0, the async processor will serialize all jobs.
	//In that case, QueueAsyncTask() will just run the asynchronous task itself immediately.
	//Returns 0 on success.
	int Initialize(int depth, bool process = false, size_t process_buffer_size = 0);
	void SetFullQueueWarning(int val) {fMe->fFullQueueWarning = val;}

	int Deinitialize();
	//Deinitialize async queue. Terminates async thread and frees buffers.
	//Requires that no tasks are left in queue. Call WaitForTasks and RetrieveQueuedTaskResult before.
	//Return 1 if it cannot deinitialize.
	
	//Returns the queue depth
	int GetQueueDepth() {return fMe->fQueueDepth;}

	//Returns the number of total async tasks currently in the queue
	int GetNumberOfAsyncTasksInQueue();

	//Can be used to perform additional initialization inside the asynchronously spawned thread.
	//Pass the initialization function as argument.
	//InitializeAsyncTask will wait for the passed initialization function in the async thread to finish.
	//This second initialization step is optional, and it is needed if some thread-local contexts must be initialized.
	//Returns the return value returned by the passed function in the async thread.
	//Must not be called whily async tasks are running, will return NULL in that case
	//In principle, the interface is the same as QueueAsyncTask below, but the function will wait for the task
	//to finish, collect the result, and return it automatically.
	int InitializeAsyncTask(void* (*initFunction)(void*), void* data, void** pRetVal = NULL);

	//Queue an asynchronous task. Provide a function pointer to the processing function and a pointer to a data object.
	//The asynchronous task shall take ownership of this data object, so the main task may no longer use it.
	//Returns 0 on success. Fails if the queue limit has already been reached.
	int QueueAsyncTask(void* (*)(void*), void*);

	//Returns the number of async tasks that are completed.
	//If the return value is greater than 1, RetrieveQueuedTaskResult() can be used to fetch the result of the first task that finished.
	int IsQueuedTaskCompleted();

	//Fetch the result of the first task that finished.
	//Returns NULL if no task has finished.
	void* RetrieveQueuedTaskResult();
	
	//Force the child process to quit. Only works with async process not thread.
	//WaitTime is the maximum time allowed for the child to quit properly in msec.
	//Returns 0 if the child exited properly, 1 if the child was killed hard, -1 on error
	int ForceChildExit(int waitTime);

	//Wait for n async tasks to finish, where n is the argument passed.
	//If n=0, the function waits for all queued tasks to finish.
	//If n is greater than the number of queued tasks, it will also wait for all tasks.
	void WaitForTasks(int);

	//The mutex can be used to serialize execution of the main task and the asynchronous task for non-thread-save parts.
	//This should be used carefully.
	//In particular, the main task must not be locked out for a long time. Return 0 on success
	int LockMutex(); //Gets the lock on the mutex
	int UnlockMutex(); //Relieves the lock on the mutex
	int TryLockMutex(); //Tries to get the lock on the mutex. Returns 0 on success if mutex could be locked. Always returns immediately.
	
	//Get size of async process shared buffer objects
	int GetBufferSize() {return(fMe->fBufferSize);}
	
	struct AliHLTAsyncProcessorBuffer
	{
		friend class AliHLTAsyncProcessor;

		void* fPtr;
		size_t fSize;

		private:
		AliHLTAsyncProcessorBuffer* fNext;
		bool fTObject;
	};
	static const size_t fgkBufferHeaderSize = sizeof(AliHLTAsyncProcessorBuffer) + (ALIHLTASYNCPROCESSOR_ALIGN - sizeof(AliHLTAsyncProcessorBuffer) % ALIHLTASYNCPROCESSOR_ALIGN) % ALIHLTASYNCPROCESSOR_ALIGN;
	
	struct AliHLTAsyncProcessorMultiBuffer
	{
		size_t fNumberOfEntries;
		AliHLTAsyncProcessorBuffer* fFirst;
	};
	static const size_t fgkMultiBufferHeaderSize = sizeof(AliHLTAsyncProcessorMultiBuffer) + (ALIHLTASYNCPROCESSOR_ALIGN - sizeof(AliHLTAsyncProcessorMultiBuffer) % ALIHLTASYNCPROCESSOR_ALIGN) % ALIHLTASYNCPROCESSOR_ALIGN;
	
	//Simple version to allocate and delete a void* ptr to a buffer for an async process only
	void* AllocateBuffer();
	void FreeBuffer(void* ptr);
	
	//More elaborate version: Allocate and delete a buffer object that contains also the size (good for passing with the functions).
	//Works for async threads and async processes. For an async process, size = 0 allocates the maximum buffer size.
	AliHLTAsyncProcessorBuffer* AllocateBuffer(size_t size);
	void FreeBuffer(AliHLTAsyncProcessorBuffer* buffer);
	
	//Serializes an object into an AliHLTAsyncProcessorBuffer and returns the pointer to it, containing the pointer to the serialized object and its size
	AliHLTAsyncProcessorBuffer* SerializeIntoBuffer(TObject* obj, AliHLTComponent* cls, AliHLTAsyncProcessorMultiBuffer* multiBuf = NULL);
	
	//Allocate and Free MultiBuffer that stores multiple objects
	AliHLTAsyncProcessorMultiBuffer* AllocateMultiBuffer();
	void FreeBuffer(AliHLTAsyncProcessorMultiBuffer* ptr);
	
	//Add a new buffer to a multibuffer. If ptr is nonzero, copy the content of ptr there.
	AliHLTAsyncProcessorBuffer* AddBuffer(AliHLTAsyncProcessorMultiBuffer* multiBuf, size_t size, void* ptr = NULL);
	AliHLTAsyncProcessorBuffer* GetEntry(AliHLTAsyncProcessorMultiBuffer* multiBuf, int num);
	
	size_t BufferSize() {return(fMe->fBufferSize);}
	void RequestPush() {fMe->fRequestPush = true;}
	bool PushRequested() {bool tmp = fMe->fRequestPush;fMe->fRequestPush = false;return tmp;}

private:
	AliHLTAsyncProcessor(const AliHLTAsyncProcessor&);
	AliHLTAsyncProcessor& operator=(const AliHLTAsyncProcessor&);
	
	//Provide additional shared buffer resources for a derived class when mode is fAsyncProcess.
	//Memory is zero-initialized.
	virtual size_t ChildSharedProcessBufferSize();
	
	size_t GetTotalSize(AliHLTAsyncProcessorMultiBuffer* multiBuf);

	struct AliHLTAsyncProcessorInput
	{
		void* (*fFunction)(void*);
		void* fData;
	};

	int GetTotalQueue();
	void AsyncThread();
	static void* AsyncThreadStartHelper(void*);
	static void* AsyncThreadStop(void*);

	static void* alignPointer(void* ptr, size_t size);
	static size_t alignSize(size_t size);
	
protected:
	struct AliHLTAsyncProcessorContent
	{
		int fQueueDepth;
		bool fAsyncThreadRunning, fAsyncThreadProcessing;
		bool fExit;
		bool fChildStopped;
		AliHLTAsyncProcessorBackend* fBackend;

		AliHLTAsyncProcessorInput* fInputQueue;
		void** fOutputQueue;
		int fInputQueueUsed, fOutputQueueUsed;
		int fWaitingForTasks;
		int fFullQueueWarning;

		void* fSynchronousOutput;	//In synchronous mode, we need one output buffer, without initialization
		
		void* fBasePtr;
		void* fBufferPtr;
		bool* fBufferUsed;
		size_t fBufferSize;
		size_t fmmapSize;
		
		int fAsyncProcess;
		void* fChildBufferSpace;
		
		bool fRequestPush;
	};
	
	int LockMutex(int i);
	int UnlockMutex(int i);
	
	AliHLTAsyncProcessorContent* fMe;	//This points to the interior of the Async Processor, possible in shared memory if the async part is an individual process

	ClassDef(AliHLTAsyncProcessor, 0);
};

#endif
