#ifndef __AliCommunicationThread__
#define __AliCommunicationThread__

#include "AliStorageClientThread.h"

#include <TThread.h>

class AliStorageClientThread;

class AliCommunicationThread
{
public:
    AliCommunicationThread(AliStorageClientThread *onlineReconstructionManager);
    ~AliCommunicationThread();
    
    void Kill();
private:
    AliStorageClientThread *fManager;
    
    static void* Dispatch(void *arg)
    {
        static_cast<AliCommunicationThread*>(arg)->CommunicationHandle();
        return nullptr;
    }
    void CommunicationHandle();
    TThread *fCommunicationThread;
    
    void SetStorageParams(int maxStorageSize,int maxOccupation,int removeEvents,int eventsInChunk);
};

#endif /* defined(__AliCommunicationThread__) */
