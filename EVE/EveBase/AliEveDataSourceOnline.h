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
#include "AliEveDataSource.h"

#include <TQObject.h>
#include <TThread.h>

class AliEveDataSourceOnline : public AliEveDataSource
{
public:
    AliEveDataSourceOnline(bool storageManager=false);
    ~AliEveDataSourceOnline();
        
    void StorageManagerOk();    // *SIGNAL*
    void StorageManagerDown();  // *SIGNAL*
    void EventServerOk();       // *SIGNAL*
    void EventServerDown();     // *SIGNAL*
    
private:
    TThread *fEventListenerThread;          // Subscriber's thread receiving events from online reco
    TThread *fStorageManagerWatcherThread;  // Thread periodically checking state of Storage Manager
    
    void GetNextEvent();            // Handler of fEventListenerThread
    void CheckStorageStatus();      // Handler of fStorageManagerWatcherThread

    void GotoEvent(Int_t event);    // Method to find proper way of delivering next event to be displayed
    void NextEvent();               // Setting new event from subscriber's thread

    void MarkCurrentEvent();        // Send request to Storage Manager to mark current event
    
    Bool_t ReceivePromptRecoParameters(Int_t runNo);
    
    AliESDEvent *fCurrentEvent[2];  // Array of currently displayed event and event being recently received from reco
    int fEventInUse;                // Specifies which element of fCurrentEvent is currently displayed
    int fWritingToEventIndex;       // Specifies to which element of fCurrentEvent subsriber's thread is writing
    bool fIsNewEventAvaliable;      // New event from reco flag
    int fFailCounter;               // How many times there were no event to be displayed
    
    Bool_t fStorageDown;            // Flag specifing if communication with Storage Manager is fine
    Bool_t fFinished;               // Will be set to true in destructor
    bool fStorageManager;           // Should Storage Manager's features be activated
    AliEveEventManager *fEventManager; // Pointer to AliEveEventManager
    
    // dispatchers for threads' handlers:
    static void* DispatchEventListener(void *arg){static_cast<AliEveDataSourceOnline*>(arg)->GetNextEvent();return nullptr;}
    static void* DispatchStorageManagerWatcher(void *arg){static_cast<AliEveDataSourceOnline*>(arg)->CheckStorageStatus();return nullptr;}
    
    ClassDef(AliEveDataSourceOnline, 0); // Interface for getting all event components in a uniform way.
};


#endif
