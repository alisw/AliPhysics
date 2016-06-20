#ifndef __AliEventsCollectorThread__
#define __AliEventsCollectorThread__

#include "AliStorageDatabase.h"
#include "AliStorageClientThread.h"

#include <TThread.h>
#include <TFile.h>

class AliStorageClientThread;

class AliEventsCollectorThread
{
public:
    AliEventsCollectorThread(AliStorageClientThread *onlineReconstructionManager);
    ~AliEventsCollectorThread();
    
    void Kill();
private:
    AliStorageClientThread *fManager;
    
    static void* Dispatch(void *arg)
    {
        static_cast<AliEventsCollectorThread*>(arg)->CollectorHandle();
        return nullptr;
    }
    void CollectorHandle();
    TThread *fCollectorThread;
    
    TFile *fCurrentFile;
    AliStorageDatabase *fDatabase;
    void CheckCurrentStorageSize();
    Long64_t GetSizeOfAllChunks();
};


#endif /* defined(__AliEventsCollectorThread__) */
