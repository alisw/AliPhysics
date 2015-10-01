/* This file is property of and copyright by the ALICE HLT Project        * 
* ALICE Experiment at CERN, All rights reserved.                         *
* See cxx source for full Copyright notice                               */

/** @file    AliHLTAsyncProcessor.cxx
@author  David Rohr (drohr@cern.ch)
*/

#include "AliHLTAsyncProcessor.h"
#include "AliHLTAsyncProcessorBackend.h"

ClassImp(AliHLTAsyncProcessor)

AliHLTAsyncProcessor::AliHLTAsyncProcessor() : AliHLTLogging(),
	fQueueDepth(0),
	fAsyncThreadRunning(false),
	fAsyncThreadProcessing(false),
	fExit(false),
	fBackend(NULL),
	fInputQueue(NULL),
	fOutputQueue(NULL),
	fInputQueueUsed(0),
	fOutputQueueUsed(0),
	fWaitingForTasks(0),
	fSynchronousOutput(NULL)
{
}

AliHLTAsyncProcessor::~AliHLTAsyncProcessor()
{
	if (fBackend) delete fBackend;
}

int AliHLTAsyncProcessor::Initialize(int depth)
{
	HLTInfo("Initializing ASYNC Processor");
	if (fQueueDepth) return(1);
	fQueueDepth = depth;
	if (fQueueDepth)
	{
		fBackend = new AliHLTAsyncProcessorBackend;
		if (fBackend == NULL || fBackend->Initialize()) return(1);
		fInputQueue = new AliHLTAsyncProcessorInput[fQueueDepth];
		fOutputQueue = new void*[fQueueDepth];
		for (int i = 2;i <= 4;i++) fBackend->LockMutex(i); //Lock thread control mutexes
		fBackend->StartThread(AsyncThreadStartHelper, this);
	}
	return(0);
}

int AliHLTAsyncProcessor::Deinitialize()
{
	HLTInfo("Deinitializing ASYNC Processor");
	if (fQueueDepth == 0) return(0);
	fBackend->LockMutex(1);
	if (GetTotalQueue())
	{
		fBackend->UnlockMutex(1);
		HLTError("Error during deinitialization of ASYNC Processor - Still tasks in queue");
		return(1);
	}
	fBackend->UnlockMutex(1);
	QueueAsyncTask(AsyncThreadStop, this);
	fBackend->StopThread();
	for (int i = 3;i <= 4;i++) fBackend->UnlockMutex(i); //Unlock remaining thread control mutexes
	delete[] fInputQueue;
	delete[] fOutputQueue;
	fAsyncThreadRunning = fAsyncThreadProcessing = false;
	delete fBackend;
	fBackend = NULL;
	fQueueDepth = 0;
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
	((AliHLTAsyncProcessor*) obj)->fExit = true;
	return(NULL);
}

void AliHLTAsyncProcessor::AsyncThread()
{
	HLTInfo("Async Thread Running");
	while (true)
	{
		fBackend->LockMutex(2); //Lock async thread if nothing to do
		while (true)
		{
			fBackend->LockMutex(1);
			if (fInputQueueUsed == 0)
			{
				fAsyncThreadRunning = false;
				fBackend->UnlockMutex(1);
				break; //No work to do, finish loops and lock this thread again
			}
			void* (*function)(void*) = fInputQueue[0].fFunction;
			void* data = fInputQueue[0].fData;
			for (int i = 0;i < fInputQueueUsed - 1;i++)
			{
				fInputQueue[i] = fInputQueue[i + 1];
			}
			fInputQueueUsed--;
			fAsyncThreadProcessing = true;
			fBackend->UnlockMutex(1);
			void* retVal = function(data);
			fBackend->LockMutex(1);
			if (fExit)
			{
				fAsyncThreadProcessing = false;
				fAsyncThreadRunning = false;
				fBackend->UnlockMutex(1);
				break;
			}
			if (fOutputQueueUsed == fQueueDepth)
			{
				fBackend->UnlockMutex(1);
				fBackend->LockMutex(4);
				fBackend->LockMutex(1);
			}
			fAsyncThreadProcessing = false;
			fOutputQueue[fOutputQueueUsed++] = retVal;
			if (fWaitingForTasks && fWaitingForTasks <= fOutputQueueUsed)
			{
				fWaitingForTasks = 0;
				fBackend->UnlockMutex(3); //Enough results in queue for WaitForTasks()
			}
			fBackend->UnlockMutex(1);
		}
		if (fExit) break;
	}
	HLTInfo("Async Thread Terminating");
}

int AliHLTAsyncProcessor::GetTotalQueue()
{
	return(fInputQueueUsed + fOutputQueueUsed + (int) fAsyncThreadProcessing);
}

int AliHLTAsyncProcessor::GetNumberOfAsyncTasksInQueue()
{
	if (fQueueDepth == 0) return(0);
	fBackend->LockMutex(1);
	int retVal = GetTotalQueue();
	fBackend->UnlockMutex(1);
	return(retVal);
}

void* AliHLTAsyncProcessor::InitializeAsyncTask(void* (*initFunction)(void*), void* data)
{
	HLTInfo("Running Initialization of ASYNC Task");
	if (GetNumberOfAsyncTasksInQueue()) return(NULL);
	QueueAsyncTask(initFunction, data);
	WaitForTasks(1);
	void* retVal = RetrieveQueuedTaskResult();
	HLTInfo("Initialization of ASYNC Task finished");
	return(retVal);
}

int AliHLTAsyncProcessor::QueueAsyncTask(void* (*processFunction)(void*), void* data)
{
	if (fQueueDepth == 0)
	{
		if (fOutputQueueUsed) return(1);
		fOutputQueueUsed = 1;
		fSynchronousOutput = processFunction(data);
		return(0);
	}
	HLTInfo("Queuing task (Queue Fill Status %d %d %d)", fInputQueueUsed, (int) fAsyncThreadProcessing, fOutputQueueUsed);
	fBackend->LockMutex(1);
	if (GetTotalQueue() == fQueueDepth)
	{
		fBackend->UnlockMutex(1);
		HLTWarning("Cannot Queue Task... Queue Full");
		return(1);
	}
	fInputQueue[fInputQueueUsed].fFunction = processFunction;
	fInputQueue[fInputQueueUsed].fData = data;
	fInputQueueUsed++;
	if (!fAsyncThreadRunning)
	{
		fAsyncThreadRunning = true;
		fBackend->UnlockMutex(2);
	}
	fBackend->UnlockMutex(1);
	return(0);
}

int AliHLTAsyncProcessor::IsQueuedTaskCompleted()
{
	HLTInfo("%d results ready for retrieval", fOutputQueueUsed);
	return(fOutputQueueUsed);
}

void* AliHLTAsyncProcessor::RetrieveQueuedTaskResult()
{
	void* retVal;
	if (fQueueDepth == 0)
	{
		fOutputQueueUsed = 0;
		retVal = fSynchronousOutput;
		fSynchronousOutput = NULL;
		return(retVal);
	}
	HLTInfo("Retrieving Queued Result");
	if (fOutputQueueUsed == 0) return(NULL);
	fBackend->LockMutex(1);
	retVal = fOutputQueue[0];
	for (int i = 0;i < fOutputQueueUsed - 1;i++)
	{
		fOutputQueue[i] = fOutputQueue[i + 1];
	}
	if (fOutputQueueUsed-- == fQueueDepth) fBackend->UnlockMutex(4); //There is no space in output queue, async thread can go on
	fBackend->UnlockMutex(1);
	return(retVal);
}

void AliHLTAsyncProcessor::WaitForTasks(int n)
{
	if (fQueueDepth == 0) return;
	fBackend->LockMutex(1);
	if (n == 0) n = fQueueDepth;
	if (n > GetTotalQueue()) n = GetTotalQueue();
	HLTInfo("Waiting for %d tasks", n);
	if (n <= fOutputQueueUsed)
	{
		fBackend->UnlockMutex(1);
		HLTInfo("%d Tasks already ready, no need to wait", n);
		return;
	}
	fWaitingForTasks = n;
	fBackend->UnlockMutex(1);
	fBackend->LockMutex(3);
	HLTInfo("Waiting for %d tasks finished", n);
}

int AliHLTAsyncProcessor::LockMutex()
{
	if (fQueueDepth == 0) return(0);
	return(fBackend->LockMutex(0));
}

int AliHLTAsyncProcessor::UnlockMutex()
{
	if (fQueueDepth == 0) return(0);
	return(fBackend->UnlockMutex(0));
}

int AliHLTAsyncProcessor::TryLockMutex()
{
	if (fQueueDepth == 0) return(0);
	return(fBackend->TryLockMutex(0));
}
