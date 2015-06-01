//
//  AliEveOnlineEventManager.h
//  xAliRoot
//
//  Created by Jeremi Niedziela on 28/05/15.
//
//

#ifndef __AliEveOnlineEventManager__
#define __AliEveOnlineEventManager__

#include "AliEveEventManager.h"
#include "AliStorageTypes.h"

class AliEveOnlineEventManager : public AliEveEventManager
{
public:
    AliEveOnlineEventManager(Int_t ev=0, bool storageManager=false);
    ~AliEveOnlineEventManager();
    
    Int_t NewEventAvailable();
    
    void StorageManagerOk();    // *SIGNAL*
    void StorageManagerDown();  // *SIGNAL*
    void EventServerOk();    // *SIGNAL*
    void EventServerDown();  // *SIGNAL*
    
private:
    void GetNextEvent();
    void CheckStorageStatus();
    void InitOCDB(int runNo);
    void GotoEvent(Int_t event);
    void NextEvent();
    void MarkCurrentEvent();
    
    TThread *fEventListenerThread;
    TThread *fStorageManagerWatcherThread;
    TMutex *fMutex;
    AliESDEvent *fCurrentEvent[2];
    TTree *fCurrentTree[2];
    int fEventInUse;
    int fWritingToEventIndex;
    bool fIsNewEventAvaliable;
    int fFailCounter;
    storageSockets fgSubSock;
    int fCurrentRun;
    
    Bool_t fOnlineMode;
    Bool_t fStorageDown;
    Bool_t fEventServerDown;
    Bool_t fFinished;
    bool fStorageManager;
    
    static void* DispatchEventListener(void *arg){static_cast<AliEveOnlineEventManager*>(arg)->GetNextEvent();return nullptr;}
    static void* DispatchStorageManagerWatcher(void *arg){static_cast<AliEveOnlineEventManager*>(arg)->CheckStorageStatus();return nullptr;}
};


#endif
