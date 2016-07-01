/* This file is property of and copyright by the ALICE HLT Project        * 
* ALICE Experiment at CERN, All rights reserved.                         *
* See cxx source for full Copyright notice                               */

/** @file    AliHLTAsyncProcessor.cxx
@author  David Rohr (drohr@cern.ch)
*/

#ifdef __MACH__
#define MAP_ANONYMOUS MAP_ANON
#endif

#include "AliHLTAsyncProcessor.h"
#include "AliHLTAsyncProcessorBackend.h"
#include "AliHLTComponent.h"
#include <sys/mman.h>

ClassImp(AliHLTAsyncProcessor)

AliHLTAsyncProcessor::AliHLTAsyncProcessor() : AliHLTLogging(), fMe(new AliHLTAsyncProcessorContent)
{
	fMe->fQueueDepth = 0;
	fMe->fAsyncThreadRunning = false;
	fMe->fAsyncThreadProcessing = false;
	fMe->fExit = false;
	fMe->fChildStopped = false;
	fMe->fBackend = NULL;
	fMe->fInputQueue = NULL;
	fMe->fOutputQueue = NULL;
	fMe->fInputQueueUsed = 0;
	fMe->fOutputQueueUsed = 0;
	fMe->fWaitingForTasks = 0;
	fMe->fFullQueueWarning = 1;
	fMe->fSynchronousOutput = NULL;
	fMe->fBasePtr = NULL;
	fMe->fBufferPtr = NULL;
	fMe->fBufferSize = 0;
	fMe->fAsyncProcess = 0;
	fMe->fChildBufferSpace = NULL;
	fMe->fmmapSize = 0;
	fMe->fRequestPush = false;
}

AliHLTAsyncProcessor::~AliHLTAsyncProcessor()
{
	if (fMe->fQueueDepth) Deinitialize();
	delete fMe;
}

size_t AliHLTAsyncProcessor::alignSize(size_t size)
{
	if (size % ALIHLTASYNCPROCESSOR_ALIGN) size += ALIHLTASYNCPROCESSOR_ALIGN - size % ALIHLTASYNCPROCESSOR_ALIGN;
	return(size);
}

void* AliHLTAsyncProcessor::alignPointer(void* ptr, size_t size)
{
	size_t tmp = (size_t) ptr;
	tmp += size;
	tmp = alignSize(tmp);
	return (void*) tmp;
}

int AliHLTAsyncProcessor::Initialize(int depth, bool process, size_t process_buffer_size)
{
	HLTInfo("Initializing ASYNC Processor");
	if (fMe->fQueueDepth) return(1);
	fMe->fQueueDepth = depth;
	if (fMe->fQueueDepth)
	{
		fMe->fAsyncProcess = process;
		size_t size = sizeof(AliHLTAsyncProcessorBackend) +
		              (sizeof(AliHLTAsyncProcessorInput) + sizeof(void*)) * fMe->fQueueDepth +
		              process_buffer_size * (fMe->fQueueDepth + 3) +
		              ChildSharedProcessBufferSize() +
		              (6 + fMe->fQueueDepth) * ALIHLTASYNCPROCESSOR_ALIGN;
		void* tmpPtr;
		if (fMe->fAsyncProcess) //promote to running async process instead of async thread
		{
			fMe->fmmapSize = sizeof(AliHLTAsyncProcessorContent) +
			                 sizeof(bool) * (fMe->fQueueDepth + 3) +
			                 2 * ALIHLTASYNCPROCESSOR_ALIGN +
			                 size;
			fMe->fBasePtr = mmap(NULL, fMe->fmmapSize, PROT_READ | PROT_WRITE, MAP_SHARED | MAP_ANONYMOUS, 0, 0);
			if (fMe->fBasePtr == NULL) return(1);
			memcpy(fMe->fBasePtr, fMe, sizeof(AliHLTAsyncProcessorContent));
			delete fMe;
			fMe = (AliHLTAsyncProcessorContent*) fMe->fBasePtr;
			fMe->fBufferSize = process_buffer_size;
			if (fMe->fBufferSize % ALIHLTASYNCPROCESSOR_ALIGN) fMe->fBufferSize += ALIHLTASYNCPROCESSOR_ALIGN - fMe->fBufferSize % ALIHLTASYNCPROCESSOR_ALIGN;
			tmpPtr = alignPointer(fMe->fBasePtr, sizeof(AliHLTAsyncProcessorContent));
			
			fMe->fBufferUsed = new (tmpPtr) bool[fMe->fQueueDepth + 3];
			memset(fMe->fBufferUsed, 0, sizeof(bool) * (fMe->fQueueDepth + 3));
			tmpPtr = alignPointer(tmpPtr, sizeof(bool) * (fMe->fQueueDepth + 3));
		}
		else
		{
			fMe->fBasePtr = malloc(size);
			if (fMe->fBasePtr == 0) return(1);
			tmpPtr = fMe->fBasePtr;
		}
		
		fMe->fBackend = new (tmpPtr) AliHLTAsyncProcessorBackend;
		if (fMe->fBackend->Initialize(fMe->fAsyncProcess)) return(1);
		tmpPtr = alignPointer(tmpPtr, sizeof(AliHLTAsyncProcessorBackend));
		
		fMe->fInputQueue = new (tmpPtr) AliHLTAsyncProcessorInput[fMe->fQueueDepth];
		tmpPtr = alignPointer(tmpPtr, sizeof(AliHLTAsyncProcessorInput) * fMe->fQueueDepth);
		
		fMe->fOutputQueue = new (tmpPtr) void*[fMe->fQueueDepth];
		tmpPtr = alignPointer(tmpPtr, sizeof(void*) * fMe->fQueueDepth);
		
		if (ChildSharedProcessBufferSize())
		{
			fMe->fChildBufferSpace = tmpPtr;
			memset(fMe->fChildBufferSpace, 0, ChildSharedProcessBufferSize());
			tmpPtr = alignPointer(tmpPtr, ChildSharedProcessBufferSize());
		}
		
		fMe->fBufferPtr = tmpPtr;
		
		for (int i = 2;i <= 4;i++) fMe->fBackend->LockMutex(i); //Lock thread control mutexes
		if (fMe->fAsyncProcess)
		{
			if (fMe->fBackend->StartProcess(AsyncThreadStartHelper, this)) return(1);
		}
		else
		{
			if (fMe->fBackend->StartThread(AsyncThreadStartHelper, this)) return(1);
		}
	}
	return(0);
}

int AliHLTAsyncProcessor::Deinitialize()
{
	HLTInfo("Deinitializing ASYNC Processor");
	if (fMe->fQueueDepth == 0) return(0);
	fMe->fBackend->LockMutex(1);
	if (GetTotalQueue())
	{
		fMe->fBackend->UnlockMutex(1);
		HLTError("Error during deinitialization of ASYNC Processor - Still tasks in queue");
		return(1);
	}
	fMe->fBackend->UnlockMutex(1);
	if (!fMe->fChildStopped) QueueAsyncTask(AsyncThreadStop, this);
	if (fMe->fAsyncProcess)
	{
		fMe->fBackend->StopProcess();
	}
	else
	{
		fMe->fBackend->StopThread();
	}
	for (int i = 3;i <= 4;i++) fMe->fBackend->UnlockMutex(i); //Unlock remaining thread control mutexes
	for (int i = 0;i < fMe->fQueueDepth;i++)
	{
		fMe->fInputQueue[i].~AliHLTAsyncProcessorInput();
		//fMe->fOutputQueue.~void*(); //Adapt if this is becoming a real struct!
	}
	fMe->fAsyncThreadRunning = fMe->fAsyncThreadProcessing = false;
	fMe->fBackend->~AliHLTAsyncProcessorBackend();;
	fMe->fBackend = NULL;
	fMe->fQueueDepth = 0;
	fMe->fBufferSize = 0;
	if (fMe->fAsyncProcess)
	{
		void* tmp = fMe;
		fMe = new AliHLTAsyncProcessorContent;
		memcpy(fMe, tmp, sizeof(AliHLTAsyncProcessorContent));
		munmap(fMe->fBasePtr, fMe->fmmapSize);
	}
	else
	{
		free(fMe->fBasePtr);
	}

	fMe->fExit = false;
	fMe->fChildStopped = false;
	fMe->fBackend = NULL;
	fMe->fInputQueue = NULL;
	fMe->fOutputQueue = NULL;
	fMe->fInputQueueUsed = 0;
	fMe->fOutputQueueUsed = 0;
	fMe->fWaitingForTasks = 0;
	fMe->fBasePtr = NULL;
	fMe->fBufferPtr = NULL;
	fMe->fBufferSize = 0;
	fMe->fAsyncProcess = 0;
	fMe->fChildBufferSpace = NULL;
	fMe->fmmapSize = 0;
	fMe->fRequestPush = false;

	HLTInfo("Deinitialization of ASYNC Processor done");
	return(0);
}

void* AliHLTAsyncProcessor::AsyncThreadStartHelper(void* obj)
{
	((AliHLTAsyncProcessor*) obj)->AsyncThread();
	pthread_exit(NULL);
	return(NULL);
}

void* AliHLTAsyncProcessor::AsyncThreadStop(void* obj)
{
	((AliHLTAsyncProcessor*) obj)->fMe->fExit = true;
	return(NULL);
}

void AliHLTAsyncProcessor::AsyncThread()
{
	HLTInfo("Async Thread Running");
	while (true)
	{
		fMe->fBackend->LockMutex(2); //Lock async thread if nothing to do
		while (true)
		{
			fMe->fBackend->LockMutex(1);
			if (fMe->fInputQueueUsed == 0)
			{
				fMe->fAsyncThreadRunning = false;
				fMe->fBackend->UnlockMutex(1);
				break; //No work to do, finish loops and lock this thread again
			}
			void* (*function)(void*) = fMe->fInputQueue[0].fFunction;
			void* data = fMe->fInputQueue[0].fData;
			for (int i = 0;i < fMe->fInputQueueUsed - 1;i++)
			{
				fMe->fInputQueue[i] = fMe->fInputQueue[i + 1];
			}
			fMe->fInputQueueUsed--;
			fMe->fAsyncThreadProcessing = true;
			fMe->fBackend->UnlockMutex(1);
			void* retVal = function(data);
			fMe->fBackend->LockMutex(1);
			if (fMe->fExit)
			{
				fMe->fAsyncThreadProcessing = false;
				fMe->fAsyncThreadRunning = false;
				fMe->fBackend->UnlockMutex(1);
				break;
			}
			if (fMe->fOutputQueueUsed == fMe->fQueueDepth)
			{
				fMe->fBackend->UnlockMutex(1);
				fMe->fBackend->LockMutex(4);
				fMe->fBackend->LockMutex(1);
			}
			fMe->fAsyncThreadProcessing = false;
			fMe->fOutputQueue[fMe->fOutputQueueUsed++] = retVal;
			if (fMe->fWaitingForTasks && fMe->fWaitingForTasks <= fMe->fOutputQueueUsed)
			{
				fMe->fWaitingForTasks = 0;
				fMe->fBackend->UnlockMutex(3); //Enough results in queue for WaitForTasks()
			}
			fMe->fBackend->UnlockMutex(1);
		}
		if (fMe->fExit) break;
	}
	HLTInfo("Async Thread Terminating");
}

int AliHLTAsyncProcessor::GetTotalQueue()
{
	return(fMe->fInputQueueUsed + fMe->fOutputQueueUsed + (int) fMe->fAsyncThreadProcessing);
}

int AliHLTAsyncProcessor::GetNumberOfAsyncTasksInQueue()
{
	if (fMe->fQueueDepth == 0) return(0);
	fMe->fBackend->LockMutex(1);
	int retVal = GetTotalQueue();
	fMe->fBackend->UnlockMutex(1);
	return(retVal);
}

int AliHLTAsyncProcessor::InitializeAsyncTask(void* (*initFunction)(void*), void* data, void** pRetVal)
{
	HLTInfo("Running Initialization of ASYNC Task");
	if (GetNumberOfAsyncTasksInQueue()) return(1);
	QueueAsyncTask(initFunction, data);
	WaitForTasks(1);
	void* retVal = RetrieveQueuedTaskResult();
	HLTInfo("Initialization of ASYNC Task finished");
	if (pRetVal) *pRetVal = retVal;
	return(0);
}

int AliHLTAsyncProcessor::QueueAsyncTask(void* (*processFunction)(void*), void* data)
{
	if (fMe->fChildStopped)
	{
		HLTError("Cannot queue new tasks after the child was stopped");
		return(1);
	}
	if (fMe->fQueueDepth == 0)
	{
		if (fMe->fOutputQueueUsed) return(1);
		fMe->fOutputQueueUsed = 1;
		fMe->fSynchronousOutput = processFunction(data);
		return(0);
	}
	HLTInfo("Queuing task (Queue Fill Status %d %d %d)", fMe->fInputQueueUsed, (int) fMe->fAsyncThreadProcessing, fMe->fOutputQueueUsed);
	fMe->fBackend->LockMutex(1);
	if (GetTotalQueue() == fMe->fQueueDepth)
	{
		fMe->fBackend->UnlockMutex(1);
		if (fMe->fFullQueueWarning) HLTWarning("Cannot Queue Task... Queue Full");
		return(1);
	}
	fMe->fInputQueue[fMe->fInputQueueUsed].fFunction = processFunction;
	fMe->fInputQueue[fMe->fInputQueueUsed].fData = data;
	fMe->fInputQueueUsed++;
	if (!fMe->fAsyncThreadRunning)
	{
		fMe->fAsyncThreadRunning = true;
		fMe->fBackend->UnlockMutex(2);
	}
	fMe->fBackend->UnlockMutex(1);
	return(0);
}

int AliHLTAsyncProcessor::IsQueuedTaskCompleted()
{
	//HLTInfo("%d results ready for retrieval", fMe->fOutputQueueUsed);
	return(fMe->fOutputQueueUsed);
}

void* AliHLTAsyncProcessor::RetrieveQueuedTaskResult()
{
	void* retVal;
	if (fMe->fQueueDepth == 0)
	{
		fMe->fOutputQueueUsed = 0;
		retVal = fMe->fSynchronousOutput;
		fMe->fSynchronousOutput = NULL;
		return(retVal);
	}
	HLTInfo("Retrieving Queued Result");
	if (fMe->fOutputQueueUsed == 0) return(NULL);
	fMe->fBackend->LockMutex(1);
	retVal = fMe->fOutputQueue[0];
	for (int i = 0;i < fMe->fOutputQueueUsed - 1;i++)
	{
		fMe->fOutputQueue[i] = fMe->fOutputQueue[i + 1];
	}
	if (fMe->fOutputQueueUsed-- == fMe->fQueueDepth) fMe->fBackend->UnlockMutex(4); //There is no space in output queue, async thread can go on
	fMe->fBackend->UnlockMutex(1);
	return(retVal);
}

void AliHLTAsyncProcessor::WaitForTasks(int n)
{
	if (fMe->fQueueDepth == 0) return;
	fMe->fBackend->LockMutex(1);
	if (n == 0) n = fMe->fQueueDepth;
	if (n > GetTotalQueue()) n = GetTotalQueue();
	HLTInfo("Waiting for %d tasks", n);
	if (n <= fMe->fOutputQueueUsed)
	{
		fMe->fBackend->UnlockMutex(1);
		HLTInfo("%d Tasks already ready, no need to wait", n);
		return;
	}
	fMe->fWaitingForTasks = n;
	fMe->fBackend->UnlockMutex(1);
	fMe->fBackend->LockMutex(3);
	HLTInfo("Waiting for %d tasks finished", n);
}

int AliHLTAsyncProcessor::LockMutex()
{
	if (fMe->fQueueDepth == 0) return(0);
	return(fMe->fBackend->LockMutex(0));
}

int AliHLTAsyncProcessor::UnlockMutex()
{
	if (fMe->fQueueDepth == 0) return(0);
	return(fMe->fBackend->UnlockMutex(0));
}

int AliHLTAsyncProcessor::TryLockMutex()
{
	if (fMe->fQueueDepth == 0) return(0);
	return(fMe->fBackend->TryLockMutex(0));
}

void* AliHLTAsyncProcessor::AllocateBuffer()
{
	if (!fMe->fAsyncProcess) return NULL;
	fMe->fBackend->LockMutex(5);
	for (int i = 0;i < fMe->fQueueDepth + 3;i++)
	{
		if (!fMe->fBufferUsed[i])
		{
			fMe->fBufferUsed[i] = true;
			fMe->fBackend->UnlockMutex(5);
			return(((char*) fMe->fBufferPtr) + i * fMe->fBufferSize);
		}
	}
	fMe->fBackend->UnlockMutex(5);
	return(NULL);
}

void AliHLTAsyncProcessor::FreeBuffer(void* ptr)
{
	if (fMe->fAsyncProcess)
	{
		for (int i = 0;i < fMe->fQueueDepth + 3;i++)
		{
			if (((char*) fMe->fBufferPtr) + i * fMe->fBufferSize == (char*) ptr)
			{
				fMe->fBufferUsed[i] = false;
				return;
			}
		}
	}
}

AliHLTAsyncProcessor::AliHLTAsyncProcessorBuffer* AliHLTAsyncProcessor::AllocateBuffer(size_t size)
{
	AliHLTAsyncProcessorBuffer* retVal;
	if (fMe->fAsyncProcess)
	{
		if (size == 0) size = fMe->fBufferSize - fgkBufferHeaderSize;
		if (size + fgkBufferHeaderSize > fMe->fBufferSize) return(NULL);
		retVal = (AliHLTAsyncProcessorBuffer*) AllocateBuffer();
	}
	else
	{
		retVal = (AliHLTAsyncProcessorBuffer*) malloc(size + fgkBufferHeaderSize);
	}
	if (retVal == NULL) return(NULL);
	retVal->fSize = size;
	retVal->fPtr = (AliHLTAsyncProcessorBuffer*) (((char*) retVal) + fgkBufferHeaderSize);
	return(retVal);
}

void AliHLTAsyncProcessor::FreeBuffer(AliHLTAsyncProcessor::AliHLTAsyncProcessorBuffer* buffer)
{
	if (fMe->fAsyncProcess)
	{
		FreeBuffer((void*) buffer);
	}
	else
	{
		free(buffer);
	}
}


AliHLTAsyncProcessor::AliHLTAsyncProcessorBuffer* AliHLTAsyncProcessor::SerializeIntoBuffer(TObject* obj, AliHLTComponent* cls, AliHLTAsyncProcessor::AliHLTAsyncProcessorMultiBuffer* multiBuf)
{
	if (multiBuf)
	{
		HLTFatal("Not yet implemented!!!");
	}
	AliHLTAsyncProcessorBuffer* retVal;
	if (fMe->fAsyncProcess)
	{
		if ((retVal = AllocateBuffer(0)) == NULL) return(NULL);
		if (cls->SerializeObject(obj, retVal->fPtr, retVal->fSize))
		{
			FreeBuffer(retVal);
			return(NULL);
		}
	}
	else
	{
		size_t size = fgkBufferHeaderSize;
		char* buffer = NULL;
		if (cls->SerializeObject(obj, (void*&) buffer, size)) return(NULL);
		retVal = (AliHLTAsyncProcessorBuffer*) buffer;
		retVal->fPtr = buffer + fgkBufferHeaderSize;
		retVal->fSize = size;
	}
	return(retVal);
}

AliHLTAsyncProcessor::AliHLTAsyncProcessorBuffer* AliHLTAsyncProcessor::AddBuffer(AliHLTAsyncProcessor::AliHLTAsyncProcessorMultiBuffer* multiBuf, size_t size, void* ptr)
{
	AliHLTAsyncProcessorBuffer* retVal;
	if (fMe->fAsyncProcess)
	{
		size_t totalSize = GetTotalSize(multiBuf);
		if (totalSize + fgkBufferHeaderSize + size > fMe->fBufferSize) return(NULL);
		retVal = (AliHLTAsyncProcessorBuffer*) (((char*) multiBuf) + totalSize);
	}
	else
	{
		retVal = (AliHLTAsyncProcessorBuffer*) malloc(size + fgkBufferHeaderSize);
		if (retVal == NULL) return(NULL);
	}
	AliHLTAsyncProcessorBuffer** b = &multiBuf->fFirst;
	for (int i = 0;i < multiBuf->fNumberOfEntries;i++) b = &((*b)->fNext);
	*b = retVal;
	retVal->fSize = size;
	retVal->fPtr = (AliHLTAsyncProcessorBuffer*) (((char*) retVal) + fgkBufferHeaderSize);
	retVal->fNext = NULL;
	if (ptr) memcpy(retVal->fPtr, ptr, size);
	multiBuf->fNumberOfEntries++;
	return(retVal);
}

AliHLTAsyncProcessor::AliHLTAsyncProcessorBuffer* AliHLTAsyncProcessor::GetEntry(AliHLTAsyncProcessor::AliHLTAsyncProcessorMultiBuffer* multiBuf, int num)
{
	if (num >= multiBuf->fNumberOfEntries) return(NULL);
	AliHLTAsyncProcessorBuffer* b = multiBuf->fFirst;
	for (int i = 0;i < multiBuf->fNumberOfEntries - 1;i++) b = b->fNext;
	return(b);
}

size_t AliHLTAsyncProcessor::GetTotalSize(AliHLTAsyncProcessor::AliHLTAsyncProcessorMultiBuffer* multiBuf)
{
	size_t totalSize = fgkMultiBufferHeaderSize;
	AliHLTAsyncProcessorBuffer* b = multiBuf->fFirst;
	for (int i = 0;i < multiBuf->fNumberOfEntries;i++)
	{
		totalSize += fgkBufferHeaderSize;
		totalSize += alignSize(b->fSize);
		b = b->fNext;
	}
	return(totalSize);
}

AliHLTAsyncProcessor::AliHLTAsyncProcessorMultiBuffer* AliHLTAsyncProcessor::AllocateMultiBuffer()
{
	AliHLTAsyncProcessorMultiBuffer* retVal;
	if (fMe->fAsyncProcess)
	{
		if (fgkMultiBufferHeaderSize > fMe->fBufferSize) return(NULL);
		retVal = (AliHLTAsyncProcessorMultiBuffer*) AllocateBuffer();
	}
	else
	{
		retVal = (AliHLTAsyncProcessorMultiBuffer*) malloc(fgkBufferHeaderSize);
	}
	if (retVal == NULL) return(NULL);
	retVal->fNumberOfEntries = 0;
	retVal->fFirst = NULL;
	return(retVal);
}

void AliHLTAsyncProcessor::FreeBuffer(AliHLTAsyncProcessor::AliHLTAsyncProcessorMultiBuffer* buffer)
{
	if (fMe->fAsyncProcess)
	{
		FreeBuffer((void*) buffer);
	}
	else
	{
		AliHLTAsyncProcessorBuffer* b = buffer->fFirst;
		for (int i = 0;i < buffer->fNumberOfEntries;i++)
		{
			AliHLTAsyncProcessorBuffer* del = b;
			b = b->fNext;
			FreeBuffer(del);
		}
		free(buffer);
	}
}

size_t AliHLTAsyncProcessor::ChildSharedProcessBufferSize()
{
	return 0;
}

int AliHLTAsyncProcessor::LockMutex(int i) {return(fMe->fBackend->LockMutex(i));}
int AliHLTAsyncProcessor::UnlockMutex(int i) {return(fMe->fBackend->UnlockMutex(i));}

int AliHLTAsyncProcessor::ForceChildExit(int waitTime)
{
	if (fMe->fChildStopped) return(0);
	if (!fMe->fAsyncProcess) return(-1);
	if (fMe->fQueueDepth == 0) return(0);
	fMe->fExit = true;
	QueueAsyncTask(AsyncThreadStop, this);
	int retVal = fMe->fBackend->KillChildProcess(waitTime);
	if (retVal == 1) HLTWarning("Async Worker Child dit not terminate in time and was killed");
	fMe->fChildStopped = true;
	fMe->fInputQueueUsed = 0;
	fMe->fAsyncThreadProcessing = false;
	fMe->fAsyncThreadRunning = false;
	if (retVal == -1) HLTError("Error: Cannot stop async child");
	return(retVal);
}
