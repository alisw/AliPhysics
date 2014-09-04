#ifndef AliStorageClientThread_H
#define AliStorageClientThread_H

#include "AliStorageDatabase.h"

#include <string>

#include <TThread.h>

class AliStorageClientThread
{
public:
	AliStorageClientThread();
	 ~AliStorageClientThread();
	void CollectData();

private:
	//status flags
	Int_t fConnectionStatus;
	Int_t fReceivingStatus;
	Int_t fSavingStatus;
	
	//communication with admin panel
    static void* Dispatch(void *arg){static_cast<AliStorageClientThread*>(arg)->CommunicationHandle();}
	void CommunicationHandle();
	TThread *fCommunicationThread;
	
	//storage file system
	void CheckCurrentStorageSize();
	void SetStorageParams(int maxStorageSize,
			      int maxOccupation,
			      int removeEvents,
			      int eventsInChunk);
	TFile *fCurrentFile;

	AliStorageDatabase *fDatabase;
	int fCurrentStorageSize;
	int fMaximumStorageSize;
	std::string fStoragePath;
	int fNumberOfEventsInFile;
	int fStorageOccupationLevel;
	int fRemoveEventsPercentage;

	Long64_t GetSizeOfAllChunks();
	
	AliStorageClientThread(const AliStorageClientThread&);
	AliStorageClientThread& operator=(const AliStorageClientThread&);
	
};

#endif
