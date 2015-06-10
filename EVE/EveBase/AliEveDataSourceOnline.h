//
//  AliEveDataSourceOnline.h
//  xAliRoot
//
//  Created by Jeremi Niedziela on 28/05/15.
//
//

#ifndef __AliEveDataSourceOnline__
#define __AliEveDataSourceOnline__

#include "AliEveEventManager.h"
#include "AliStorageTypes.h"

class AliEveDataSourceOnline : public AliEveDataSource
{
public:
    AliEveDataSourceOnline(bool storageManager=false);
    ~AliEveDataSourceOnline();
        
    void StorageManagerOk();    // *SIGNAL*
    void StorageManagerDown();  // *SIGNAL*
    void EventServerOk();    // *SIGNAL*
    void EventServerDown();  // *SIGNAL*
    
private:
    void GetNextEvent();
    void CheckStorageStatus();
    void Init(){};// to implement
    void InitOCDB(int runNo);
    void GotoEvent(Int_t event);
    void NextEvent();
    void MarkCurrentEvent();
    
    TThread *fEventListenerThread;
    TThread *fStorageManagerWatcherThread;
    AliESDEvent *fCurrentEvent[2];
    TTree *fCurrentTree[2];
    int fEventInUse;
    int fWritingToEventIndex;
    bool fIsNewEventAvaliable;
    int fFailCounter;
    int fCurrentRun;
    
    Bool_t fOnlineMode;
    Bool_t fStorageDown;
    Bool_t fEventServerDown;
    Bool_t fFinished;
    bool fStorageManager;
    
    static void* DispatchEventListener(void *arg){static_cast<AliEveDataSourceOnline*>(arg)->GetNextEvent();return nullptr;}
    static void* DispatchStorageManagerWatcher(void *arg){static_cast<AliEveDataSourceOnline*>(arg)->CheckStorageStatus();return nullptr;}
    
    AliEveDataSourceOnline(const AliEveDataSourceOnline&);
    AliEveDataSourceOnline& operator=(const AliEveDataSourceOnline&);
};


#endif
