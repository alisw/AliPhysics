//
//  AliEveDataSourceOnline.cpp
//  xAliRoot
//
//  Created by Jeremi Niedziela on 28/05/15.
//
//

#include "AliEveDataSourceOnline.h"
#include "AliEveConfigManager.h"
#include "AliEveInit.h"

#include "AliZMQManager.h"
#include "AliGRPPreprocessor.h"
#include "AliEveInit.h"
#include "AliCDBManager.h"

#include <TEnv.h>
#include <TInterpreter.h>

#include <iostream>

using namespace std;

ClassImp(AliEveDataSourceOnline)

AliEveDataSourceOnline::AliEveDataSourceOnline(bool storageManager) :
AliEveDataSource(storageManager),
fEventListenerThread(0),
fStorageManagerWatcherThread(0),
fEventInUse(1),
fWritingToEventIndex(0),
fIsNewEventAvaliable(false),
fFailCounter(0),
fStorageDown(false),
fHasEventFromStorageManager(false),
fFinished(false),
fStorageManager(storageManager),
fEventManager(0)
{
    fEventManager = AliEveEventManager::Instance();
    
    StorageManagerDown(); // turn SM off by default
    
    // start threads:
    fEventListenerThread = new TThread("fEventListenerThread",DispatchEventListener,(void*)this);
    fEventListenerThread->Run();
    
    if(fStorageManager){
        fStorageManagerWatcherThread = new TThread("fStorageManagerWatcherThread",DispatchStorageManagerWatcher,(void*)this);
        fStorageManagerWatcherThread->Run();
    }
}

AliEveDataSourceOnline::~AliEveDataSourceOnline()
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
}

void AliEveDataSourceOnline::GetNextEvent()
{
    AliZMQManager *eventManager = AliZMQManager::GetInstance();
    eventManager->CreateSocket(EVENTS_SERVER_SUB);
    
    fCurrentEvent[0]=0;
    fCurrentEvent[1]=0;
    
    AliESDEvent *tmpEvent;
    
    cout<<"Starting subscriber's loop"<<endl;
    int receiveStatus;
    
    while(1)
    {
        cout<<"Waiting for event from online reconstruction..."<<endl;
        receiveStatus = eventManager->Get(tmpEvent,EVENTS_SERVER_SUB);
        
        if(receiveStatus == 0){ // timeout
            continue;
        }
        else if(receiveStatus == -1){ // error, socket closed
            break;
        }
        else if(receiveStatus == 1)
        {
            if(tmpEvent)
            {
                cout<<"Received event of run "<<tmpEvent->GetRunNumber()<<endl;
                if(tmpEvent->GetRunNumber()>=0 && tmpEvent->GetNumberOfTracks()>0)
                {
#if ROOT_VERSION_CODE < ROOT_VERSION(5,99,0)
                    gCINTMutex->Lock();
#endif
                    TEnv settings;
                    AliEveInit::GetConfig(&settings);
                    const char *requiredTriggerClass = settings.GetValue("trigger.class.filter", "");
                    
                    bool matchesTriggerClass = false;
                    
                    
                    if(strcmp(requiredTriggerClass,"none")!=0)
                    {
                        // check if required trigger class fired
                        ULong64_t mask = 1;
                        ULong64_t triggerMask = tmpEvent->GetTriggerMask();
                        ULong64_t triggerMaskNext50 = tmpEvent->GetTriggerMaskNext50();
                        
                        mask=1;
                        
                        for(int i=0;i<50;i++)
                        {
                            if(mask & triggerMask)
                            {
                                if(strcmp(tmpEvent->GetESDRun()->GetTriggerClass(i),requiredTriggerClass)==0)
                                {
                                    cout<<"\n\n\nMatching trigger class\n\n\n"<<endl;
                                    matchesTriggerClass=true;
                                }
                            }
                            if(mask & triggerMaskNext50)
                                {
                                if(strcmp(tmpEvent->GetESDRun()->GetTriggerClass(i+50),requiredTriggerClass)==0)
                                {
                                    cout<<"\n\n\nMatching trigger class\n\n\n"<<endl;
                                    matchesTriggerClass=true;
                                }
                            }
                            mask = mask<<1;
                        }
                    //
                    }
                    else
                    {
                        matchesTriggerClass = true;
                    }
                    
                    if(matchesTriggerClass)
                    {
                        if(fEventInUse == 0){fWritingToEventIndex = 1;}
                        else if(fEventInUse == 1){fWritingToEventIndex = 0;}
                        cout<<"Event number in file:"<<tmpEvent->GetEventNumberInFile()<<endl;
                        if(fCurrentEvent[fWritingToEventIndex])
                        {
                            delete fCurrentEvent[fWritingToEventIndex];
                            fCurrentEvent[fWritingToEventIndex]=0;
                        }
                        fCurrentEvent[fWritingToEventIndex] = tmpEvent;
                        
                        
                        fIsNewEventAvaliable = true;
                    }
                    else
                    {
                        cout<<"Trigger classes don't match the filter!"<<endl;
                    }
#if ROOT_VERSION_CODE < ROOT_VERSION(5,99,0)
                    gCINTMutex->UnLock();
#endif
                }
                else
                {
                    if(tmpEvent->GetRunNumber()<0)
                    {
                        cout<<"Incorrect run number!"<<endl;
                    }
                    if(tmpEvent->GetNumberOfTracks()<=0)
                    {
                        cout<<"No tracks!"<<endl;
                    }
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

void AliEveDataSourceOnline::CheckStorageStatus()
{
    if(!fStorageManager){return;}
    
    AliEveConfigManager *configManager = AliEveConfigManager::Instance();
    configManager->ConnectEventManagerSignals();
    
    AliZMQManager *eventManager = AliZMQManager::GetInstance();
    storageSockets socket = CLIENT_COMMUNICATION_REQ;
    eventManager->CreateSocket(socket);
    
    struct clientRequestStruct *request = new struct clientRequestStruct;
    request->messageType = REQUEST_CONNECTION;
    
    long response;
    int receiveStatus = 0;
    int sendStatus = 0;
    
    StorageManagerDown();// assume that storage manager is down
    
    while (1)
    {
        sendStatus = eventManager->Send(request,socket);
        if(sendStatus == 0){ // timeout
            continue;
        }
        else if(sendStatus == -1){ // error, socket closed
            break;
        }
        
        cout<<"EVENT DISPLAY -- message sent to SM"<<endl;
        
        receiveStatus = eventManager->Get(&response,socket); // try to reveive response
        if(receiveStatus == 0) // timeout
        {
            cout<<"EVENT DISPLAY -- failed to receive message from SM, timeout"<<endl;
            eventManager->RecreateSocket(socket); // destroy and open socket again, to be able to send message (currently, as receive failed, socket is still in RECV state
            
            if(!fStorageDown) // if requires change
            {
                cout<<"EVENT DISPLAY -- storage DOWN"<<endl;
                StorageManagerDown();
            }
        }
        else if(receiveStatus == -1){ // error, socket closed
            cout<<"EVENT DISPLAY -- failed to receive message from SM, error"<<endl;
            break;
        }
        else if(fStorageDown) // if success and requires change
        {
            cout<<"EVENT DISPLAY -- storage OK"<<endl;
            StorageManagerOk();
        }
        sleep(1);
        sendStatus=0;
    }
    
    AliEveEventManager *manager = fEventManager;
    manager->Disconnect("StorageManagerOk");
    manager->Disconnect("StorageManagerDown");
}

void AliEveDataSourceOnline::GotoEvent(Int_t event)
{
    cout<<"Go to event:"<<event<<endl;
    
    if (fStorageDown && -1 == event)
    {
        NextEvent();
        return;
    }
    
    if (fCurrentData.fESD)
    {
        // create new server request:
        struct serverRequestStruct *requestMessage = new struct serverRequestStruct;
        
        // set request type:
        if (event == -1)      {requestMessage->messageType = REQUEST_GET_LAST_EVENT;}
        else  if (event == 0) {requestMessage->messageType = REQUEST_GET_FIRST_EVENT;}
        else  if (event == 1) {requestMessage->messageType = REQUEST_GET_PREV_EVENT;}
        else  if (event == 2) {requestMessage->messageType = REQUEST_GET_NEXT_EVENT;}
        
        // set event struct:
        requestMessage->eventsRunNumber = fCurrentData.fESD->GetRunNumber();
        requestMessage->eventsEventNumber = fCurrentData.fESD->GetEventNumberInFile();
        
        // create event manager:
        AliZMQManager *eventManager = AliZMQManager::GetInstance();
        AliESDEvent *resultEvent = NULL;
        
#if ROOT_VERSION_CODE < ROOT_VERSION(5,99,0)
        gCINTMutex->Lock();
#endif
        
        // send request and receive event:
        eventManager->Send(requestMessage,SERVER_COMMUNICATION_REQ);
        eventManager->Get(resultEvent,SERVER_COMMUNICATION_REQ);
        
        if(resultEvent)
        {
            //            DestroyElements();
            fCurrentData.fESD = resultEvent;
        }
        else
        {
            if(event==-1){cout<<"\n\nWARNING -- No last event is avaliable.\n\n"<<endl;}
            if(event==0){cout<<"\n\nWARNING -- No first event is avaliable.\n\n"<<endl;}
            if(event==1){cout<<"\n\nWARNING -- No previous event is avaliable.\n\n"<<endl;}
            if(event==2){cout<<"\n\nWARNING -- No next event is avaliable.\n\n"<<endl;}
        }
        
#if ROOT_VERSION_CODE < ROOT_VERSION(5,99,0)
        gCINTMutex->UnLock();
#endif
    }
    else
    {
        cout<<"\n\nWARNING -- No event has been already loaded. Loading the most recent event...\n\n"<<endl;
        return;
        struct serverRequestStruct *requestMessage = new struct serverRequestStruct;
        requestMessage->messageType = REQUEST_GET_LAST_EVENT;
        AliZMQManager *eventManager = AliZMQManager::GetInstance();
        AliESDEvent *resultEvent = NULL;
#if ROOT_VERSION_CODE < ROOT_VERSION(5,99,0)
        gCINTMutex->Lock();
#endif

        eventManager->Send(requestMessage,SERVER_COMMUNICATION_REQ);
        eventManager->Get(resultEvent,SERVER_COMMUNICATION_REQ);
        
        if(resultEvent)
        {
            //            DestroyElements();
            fCurrentData.fESD = resultEvent;
        }
        else{cout<<"\n\nWARNING -- The most recent event is not avaliable.\n\n"<<endl;}
#if ROOT_VERSION_CODE < ROOT_VERSION(5,99,0)
        gCINTMutex->UnLock();
#endif
    }
}

void AliEveDataSourceOnline::SetEventFromStorageManager(AliESDEvent *event)
{
    fCurrentData.fESD = event;
    fHasEventFromStorageManager = true;
}

void AliEveDataSourceOnline::NextEvent()
{
#if ROOT_VERSION_CODE < ROOT_VERSION(5,99,0)
    gCINTMutex->Lock();
#endif
    if(fHasEventFromStorageManager)
    {
        fEventManager->DestroyTransients();
        fEventManager->DestroyElements();
        
        if(fCurrentData.fESD->GetRunNumber() != fEventManager->GetCurrentRun()){
            fEventManager->ResetMagneticField();
            fEventManager->SetCurrentRun(fCurrentData.fESD->GetRunNumber());
        }
        
        fEventManager->SetHasEvent(true);
        fEventManager->AfterNewEventLoaded();
        fEventManager->NewEventLoaded();
        
        fHasEventFromStorageManager = false;
    }
    else if(fIsNewEventAvaliable)
    {
        fEventManager->DestroyTransients();
        
        if(fWritingToEventIndex == 0) fEventInUse = 0;
        else if(fWritingToEventIndex == 1) fEventInUse = 1;
        
        if(fCurrentEvent[fEventInUse])
        {
            if(fCurrentEvent[fEventInUse]->GetRunNumber() >= 0)
            {
                fFailCounter=0;
                cout<<"================ setting event to "<<fCurrentEvent[fEventInUse]->GetEventNumberInFile()<<"================"<<endl;
                
                StorageManagerDown(); // block SM while event is being loaded
                fEventManager->DestroyElements();
                if(fCurrentEvent[fEventInUse]->GetRunNumber() != fEventManager->GetCurrentRun()){
                    fEventManager->ResetMagneticField();
                    fEventManager->SetCurrentRun(fCurrentEvent[fEventInUse]->GetRunNumber());
                }
                fCurrentData.fESD = fCurrentEvent[fEventInUse];
                
                fEventManager->SetHasEvent(true);
                fEventManager->AfterNewEventLoaded();
                
                if (fEventManager->GetAutoLoad()) {
                    fEventManager->StartAutoLoadTimer();
                }
                fEventManager->NewEventLoaded();
            }
        }
        fIsNewEventAvaliable = false;
    }
    else
    {
        cout<<"No new event is avaliable."<<endl;
        fFailCounter++;
        fEventManager->NoEventLoaded();
    }
    if(fFailCounter==5)
    {
        fFailCounter=0;
    }
#if ROOT_VERSION_CODE < ROOT_VERSION(5,99,0)
    gCINTMutex->UnLock();
#endif
    gSystem->ProcessEvents();
}


void AliEveDataSourceOnline::MarkCurrentEvent()
{
    if(!fStorageManager){return;}
    
    struct serverRequestStruct *requestMessage = new struct serverRequestStruct;
    requestMessage->eventsRunNumber = fCurrentData.fESD->GetRunNumber();
    requestMessage->eventsEventNumber = fCurrentData.fESD->GetEventNumberInFile();
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

void AliEveDataSourceOnline::StorageManagerOk()
{
    //    Emit("StorageManagerOk()");
}
void AliEveDataSourceOnline::StorageManagerDown()
{
    //    Emit("StorageManagerDown()");
}

Bool_t AliEveDataSourceOnline::ReceivePromptRecoParameters(Int_t runNo)
{
  TString localGRPstorage = "local://OCDB";
  cout<<"Loading OCDB for new run:"<<runNo<<" in online mode."<<endl;
  TEnv settings;
  AliEveInit::GetConfig(&settings);

  // Retrieve GRP entry for given run from aldaqdb.
  TString dbHost = settings.GetValue("logbook.host", "");
  Int_t   dbPort =  settings.GetValue("logbook.port", 0);
  TString dbName =  settings.GetValue("logbook.db", "");
  TString user =  settings.GetValue("logbook.user", "");
  TString password = settings.GetValue("logbook.pass", "");

  gSystem->cd(localGRPstorage.Data());
  gSystem->Exec("rm -fr GRP/");
  cout<<"CDB path for GRP:"<<localGRPstorage<<endl;

  TString gdc;

  Int_t ret=AliGRPPreprocessor::ReceivePromptRecoParameters(runNo, dbHost.Data(),
                                                            dbPort, dbName.Data(),
                                                            user.Data(), password.Data(),
                                                            Form("%s",localGRPstorage.Data()),
                                                            gdc);

  if(ret>0) Info("RetrieveGRP","Last run of the same type is: %d",ret);
  else if(ret==0) Warning("RetrieveGRP","No previous run of the same type found");
  else if(ret<0) Error("Retrieve","Error code while retrieving GRP parameters returned: %d",ret);

  AliCDBManager* cdb = AliCDBManager::Instance();
  cdb->SetDefaultStorage(settings.GetValue("cdb.defaultStorage",Form("local://%s/../src/OCDB",gSystem->Getenv("ALICE_ROOT"))));
  cdb->SetSpecificStorage("GRP/GRP/Data",localGRPstorage.Data());
  cdb->SetRun(runNo);
  cdb->Print();
  return kTRUE;    
}
