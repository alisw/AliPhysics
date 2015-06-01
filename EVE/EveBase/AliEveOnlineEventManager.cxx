//
//  AliEveOnlineEventManager.cpp
//  xAliRoot
//
//  Created by Jeremi Niedziela on 28/05/15.
//
//

#include "AliEveOnlineEventManager.h"
#include "AliEveConfigManager.h"

#include "AliCDBManager.h"

#include "AliZMQManager.h"
#include "AliOnlineReconstructionUtil.h"
#include "AliGRPPreprocessor.h"

#include <TEnv.h>

#include <iostream>

using namespace std;

AliEveOnlineEventManager::AliEveOnlineEventManager(Int_t ev, bool storageManager) :
    AliEveEventManager("online",ev,true),
    fMutex(new TMutex()),
    fgSubSock(EVENTS_SERVER_SUB),
    fCurrentRun(-1),
    fEventInUse(1),
    fWritingToEventIndex(0),
    fIsNewEventAvaliable(false),
    fFailCounter(0),
    fOnlineMode(kFALSE),
    fStorageDown(false),
    fEventServerDown(false),
    fFinished(false),
    fStorageManager(storageManager)
{
    cout<<"\n\n\nAliEveOnlineEventManager constructor called!!!\n\n\n"<<endl;
    
    InitOCDB(-1);
    fIsOpen = kFALSE;
    
    StorageManagerDown(); // turn SM off by default
    EventServerDown();
    
    fEventListenerThread = new TThread("fEventListenerThread",DispatchEventListener,(void*)this);
    fEventListenerThread->Run();
    
    if(fStorageManager)
    {
        fStorageManagerWatcherThread = new TThread("fStorageManagerWatcherThread",DispatchStorageManagerWatcher,(void*)this);
        fStorageManagerWatcherThread->Run();
    }
    AliEveEventManager::SetMaster(this);
}

AliEveOnlineEventManager::~AliEveOnlineEventManager()
{
    fFinished = true;
    if(fEventListenerThread)
    {
        fEventListenerThread->Join();
        fEventListenerThread->Kill();
        delete fEventListenerThread;
        cout<<"listener thread killed and deleted"<<endl;
    }
    if(fStorageManagerWatcherThread)
    {
        fStorageManagerWatcherThread->Join();
        fStorageManagerWatcherThread->Kill();
        delete fStorageManagerWatcherThread;
        cout<<"storage watcher thread killed and deleted"<<endl;
    }
    if(fMutex){delete fMutex;}
}

void AliEveOnlineEventManager::GetNextEvent()
{
    cout<<"\n\nGet next event called\n\n"<<endl;
    
    AliZMQManager *eventManager = AliZMQManager::GetInstance();
    eventManager->CreateSocket(EVENTS_SERVER_SUB);
    
    fCurrentEvent[0]=0;
    fCurrentEvent[1]=0;
    
    AliESDEvent *tmpEvent;
    
    cout<<"Starting subscriber's loop"<<endl;
    while(!fFinished)
    {
        cout<<"Waiting for event from online reconstruction...";
        if(eventManager->Get(tmpEvent,EVENTS_SERVER_SUB))
        {
            if(tmpEvent)
            {
                cout<<"received. ("<<tmpEvent->GetRunNumber();
                if(tmpEvent->GetRunNumber()>=0)
                {
                    fMutex->Lock();
                    if(fEventInUse == 0){fWritingToEventIndex = 1;}
                    else if(fEventInUse == 1){fWritingToEventIndex = 0;}
                    cout<<","<<tmpEvent->GetEventNumberInFile()<<")"<<endl;
                    if(fCurrentEvent[fWritingToEventIndex])
                    {
                        cout<<"EVENT DISPLAY -- deleting old event..."<<fCurrentEvent[fWritingToEventIndex]<<"...";
                        delete fCurrentEvent[fWritingToEventIndex];
                        fCurrentEvent[fWritingToEventIndex]=0;
                        cout<<"deleted"<<endl;
                    }
                    fCurrentEvent[fWritingToEventIndex] = tmpEvent;
                    fIsNewEventAvaliable = true;
                    NewEventLoaded();
                    fMutex->UnLock();
                }
            }
            else
            {
                cout<<"received empty event."<<endl;
                sleep(1);
            }
        }
        else
        {
            cout<<"couldn't receive.";
            sleep(1);
        }
    }
}

void AliEveOnlineEventManager::CheckStorageStatus()
{
    if(!fStorageManager){return;}
    
    AliEveConfigManager *configManager = AliEveConfigManager::GetMaster();
    configManager->ConnectEventManagerSignals();
    
    AliZMQManager *eventManager = AliZMQManager::GetInstance();
    storageSockets socket = CLIENT_COMMUNICATION_REQ;
    eventManager->CreateSocket(socket);
    
    struct clientRequestStruct *request = new struct clientRequestStruct;
    request->messageType = REQUEST_CONNECTION;
    
    long response;
    bool receiveStatus = false;
    bool sendStatus = false;
    
    StorageManagerDown();// assume that storage manager is down
    
    while (!fFinished)
    {
        while(sendStatus==false)// send request message until success
        {
            sendStatus = eventManager->Send(request,socket);
            if(sendStatus==false)
            {
                //eventManager->RecreateSocket(socket);
                sleep(1);
            }
        }
        cout<<"EVENT DISPLAY -- message sent to SM"<<endl;
        
        receiveStatus = eventManager->Get(&response,socket); // try to reveive response
        if(receiveStatus == false) // if failed (or timeouted)
        {
            cout<<"EVENT DISPLAY -- failed to receive message from SM"<<endl;
            eventManager->RecreateSocket(socket); // destroy and open socket again, to be able to send message (currently, as receive failed, socket is still in RECV state
            
            if(!fStorageDown) // if requires change
            {
                cout<<"EVENT DISPLAY -- storage DOWN"<<endl;
                StorageManagerDown();
            }
        }
        else if(fStorageDown) // if success and requires change
        {
            cout<<"EVENT DISPLAY -- storage OK"<<endl;
            StorageManagerOk();
        }
        sleep(1);
        sendStatus=false;
    }
    
    AliEveEventManager *manager = AliEveEventManager::GetMaster();
    manager->Disconnect("StorageManagerOk");
    manager->Disconnect("StorageManagerDown");
}

void AliEveOnlineEventManager::InitOCDB(int runNo)
{
    TString cdbPath = Form("local://%s/ed_ocdb_objects/",gSystem->Getenv("HOME"));
    AliCDBManager* cdb = AliCDBManager::Instance();
    
    if(runNo != fCurrentRun)
    {
        cout<<"Loading OCDB for new run:"<<runNo<<" in online mode."<<endl;
        TEnv settings;
        settings.ReadFile(AliOnlineReconstructionUtil::GetPathToServerConf(), kEnvUser);
        fCurrentRun = runNo;
        cout<<"config read"<<endl;
        
        // Retrieve GRP entry for given run from aldaqdb.
        TString dbHost = settings.GetValue("logbook.host", DEFAULT_LOGBOOK_HOST);
        Int_t   dbPort =  settings.GetValue("logbook.port", DEFAULT_LOGBOOK_PORT);
        TString dbName =  settings.GetValue("logbook.db", DEFAULT_LOGBOOK_DB);
        TString user =  settings.GetValue("logbook.user", DEFAULT_LOGBOOK_USER);
        TString password = settings.GetValue("logbook.pass", DEFAULT_LOGBOOK_PASS);
        
        gSystem->cd(cdbPath.Data());
        gSystem->Exec("rm -fr GRP/");
        cout<<"CDB path for GRP:"<<cdbPath<<endl;
        
        TString gdc;
        
        Int_t ret=AliGRPPreprocessor::ReceivePromptRecoParameters(fCurrentRun, dbHost.Data(),
                                                                  dbPort, dbName.Data(),
                                                                  user.Data(), password.Data(),
                                                                  Form("%s",cdbPath.Data()),
                                                                  gdc);
        
        if(ret>0) Info("RetrieveGRP","Last run of the same type is: %d",ret);
        else if(ret==0) Warning("RetrieveGRP","No previous run of the same type found");
        else if(ret<0) Error("Retrieve","Error code while retrieving GRP parameters returned: %d",ret);
        
        
        cdb->SetDefaultStorage(settings.GetValue("cdb.defaultStorage", DEFAULT_CDB_STORAGE));
        cdb->SetSpecificStorage("GRP/GRP/Data",cdbPath.Data());
        cdb->SetRun(fCurrentRun);
        cdb->Print();
    }
    
    
    AliEveEventManager::InitOCDB(runNo);
}

void AliEveOnlineEventManager::GotoEvent(Int_t event)
{
    cout<<"Go to event:"<<event<<endl;
    
    static const TEveException kEH("AliEveEventManager::GotoEvent ");
    
    if (fAutoLoadTimerRunning){throw (kEH + "Event auto-load timer is running.");}
    else if (!fIsOpen && !fOnlineMode){throw (kEH + "Event-files not opened but ED is in offline mode.");}
    
    if (fStorageDown && -1 == event)
    {
        NextEvent();
        return;
    }
    
    if (fESD)
    {
        // create new server request:
        struct serverRequestStruct *requestMessage = new struct serverRequestStruct;
        
        // set request type:
        if (event == -1)      {requestMessage->messageType = REQUEST_GET_LAST_EVENT;}
        else  if (event == 0) {requestMessage->messageType = REQUEST_GET_FIRST_EVENT;}
        else  if (event == 1) {requestMessage->messageType = REQUEST_GET_PREV_EVENT;}
        else  if (event == 2) {requestMessage->messageType = REQUEST_GET_NEXT_EVENT;}
        
        // set event struct:
        requestMessage->eventsRunNumber = fESD->GetRunNumber();
        requestMessage->eventsEventNumber = fESD->GetEventNumberInFile();
        
        // create event manager:
        AliZMQManager *eventManager = AliZMQManager::GetInstance();
        AliESDEvent *resultEvent = NULL;
        
        fMutex->Lock();
        
        // send request and receive event:
        eventManager->Send(requestMessage,SERVER_COMMUNICATION_REQ);
        eventManager->Get(resultEvent,SERVER_COMMUNICATION_REQ);
        
        if(resultEvent)
        {
            DestroyElements();
            InitOCDB(resultEvent->GetRunNumber());
            SetEvent(0,0,resultEvent,0);
        }
        else
        {
            if(event==-1){cout<<"\n\nWARNING -- No last event is avaliable.\n\n"<<endl;}
            if(event==0){cout<<"\n\nWARNING -- No first event is avaliable.\n\n"<<endl;}
            if(event==1){cout<<"\n\nWARNING -- No previous event is avaliable.\n\n"<<endl;}
            if(event==2){cout<<"\n\nWARNING -- No next event is avaliable.\n\n"<<endl;}
        }
        
        fMutex->UnLock();
    }
    else
    {
        cout<<"\n\nWARNING -- No event has been already loaded. Loading the most recent event...\n\n"<<endl;
        
        struct serverRequestStruct *requestMessage = new struct serverRequestStruct;
        requestMessage->messageType = REQUEST_GET_LAST_EVENT;
        
        AliZMQManager *eventManager = AliZMQManager::GetInstance();
        AliESDEvent *resultEvent = NULL;
        
        fMutex->Lock();
        eventManager->Send(requestMessage,SERVER_COMMUNICATION_REQ);
        eventManager->Get(resultEvent,SERVER_COMMUNICATION_REQ);
        
        if(resultEvent)
        {
            DestroyElements();
            InitOCDB(resultEvent->GetRunNumber());
            SetEvent(0,0,resultEvent,0);
        }
        else{cout<<"\n\nWARNING -- The most recent event is not avaliable.\n\n"<<endl;}
        fMutex->UnLock();
    }
    //
//    AliEveEventManager::GotoEvent(event);
    //
}

void AliEveOnlineEventManager::NextEvent()
{
    static const TEveException kEH("AliEveEventManager::NextEvent ");
    
    if (fAutoLoadTimerRunning){throw (kEH + "Event auto-load timer is running.");}
    
    fMutex->Lock();
    if(fIsNewEventAvaliable)
    {
        if(fWritingToEventIndex == 0) fEventInUse = 0;
        else if(fWritingToEventIndex == 1) fEventInUse = 1;
        
        if(fCurrentEvent[fEventInUse])
        {
            if(fCurrentEvent[fEventInUse]->GetRunNumber() >= 0)
            {
                fFailCounter=0;
                EventServerOk();
                printf("======================= setting event to %d\n", fCurrentEvent[fEventInUse]->GetEventNumberInFile());
                
                StorageManagerDown(); // block SM while event is being loaded
                DestroyElements();
                InitOCDB(fCurrentEvent[fEventInUse]->GetRunNumber());
                SetEvent(0,0,fCurrentEvent[fEventInUse],0);
            }
        }
        fIsNewEventAvaliable = false;
    }
    else
    {
        cout<<"No new event is avaliable."<<endl;
        fFailCounter++;
        NoEventLoaded();
    }
    if(fFailCounter==5)
    {
        cout<<"fail counter = 5"<<endl;
        EventServerDown();
        fFailCounter=0;
    }
    fMutex->UnLock();
    
    gSystem->ProcessEvents();
}


void AliEveOnlineEventManager::MarkCurrentEvent()
{
    if(!fStorageManager){return;}
    
    struct serverRequestStruct *requestMessage = new struct serverRequestStruct;
    requestMessage->eventsRunNumber = fESD->GetRunNumber();
    requestMessage->eventsEventNumber = fESD->GetEventNumberInFile();
    requestMessage->messageType = REQUEST_MARK_EVENT;
    
    AliZMQManager *eventManager = AliZMQManager::GetInstance();
    
    eventManager->Send(requestMessage,SERVER_COMMUNICATION_REQ);
    bool response;
    eventManager->Get(&response,SERVER_COMMUNICATION_REQ);
    
    if(response)
    {
        //fStatusLabel->SetText("Event marked");
        cout<<"ADMIN PANEL -- Event marked succesfully"<<endl;
    }
    else
    {
        //fStatusLabel->SetText("Couldn't mark this event");
        cout<<"ADMIN PANEL -- Could not matk event"<<endl;
    }
    if(requestMessage){delete requestMessage;}
}

Int_t AliEveOnlineEventManager::NewEventAvailable()
{
    if (fIsNewEventAvaliable){return 1;}
    else{return 0;}
}

void AliEveOnlineEventManager::StorageManagerOk()
{
    Emit("StorageManagerOk()");
}
void AliEveOnlineEventManager::StorageManagerDown()
{
    Emit("StorageManagerDown()");
}

void AliEveOnlineEventManager::EventServerOk()
{
    Emit("EventServerOk()");
}
void AliEveOnlineEventManager::EventServerDown()
{
    Emit("EventServerDown()");
}



