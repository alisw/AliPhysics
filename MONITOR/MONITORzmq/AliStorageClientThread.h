#ifndef AliStorageClientThread_H
#define AliStorageClientThread_H

#include "AliStorageTypes.h"

#include "AliDIMListenerThread.h"
#include "AliEventsCollectorThread.h"
#include "AliCommunicationThread.h"

#include <string>

class AliCommunicationThread;
class AliEventsCollectorThread;

class AliStorageClientThread
{
    friend class AliEventsCollectorThread;
    friend class AliCommunicationThread;
    
public:
	AliStorageClientThread();
	 ~AliStorageClientThread();
private:
    AliDIMListenerThread *fDIMListenerThread;
    AliEventsCollectorThread *fEventsCollectorThread;
    AliCommunicationThread *fCommunicationThread;
    
	AliStorageClientThread(const AliStorageClientThread&);
	AliStorageClientThread& operator=(const AliStorageClientThread&);
    
protected:
    // status flags
    Int_t fConnectionStatus;
    Int_t fReceivingStatus;
    Int_t fSavingStatus;

    // storage parameters
    int fCurrentStorageSize;
    int fMaximumStorageSize;
    std::string fStoragePath;
    int fNumberOfEventsInFile;
    int fStorageOccupationLevel;
    int fRemoveEventsPercentage;
};

#endif