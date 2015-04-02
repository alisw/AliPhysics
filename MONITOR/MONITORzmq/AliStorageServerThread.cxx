#include "AliStorageServerThread.h"
#include "AliStorageTypes.h"
#include "AliESDEvent.h"

#include <iostream>
#include <fstream>

#include <TFile.h>
#include <TThread.h>

using namespace std;

AliStorageServerThread::AliStorageServerThread() :
fDatabase(0),
fStoragePath("")
{
    fDatabase = new AliStorageDatabase();
    
    TThread::Lock();
    //load parameters from config file
    ifstream configFile (GetConfigFilePath());
    if (configFile.is_open())
    {
        string line;
        int from,to;
        while(configFile.good())
        {
            getline(configFile,line);
            from = line.find("\"")+1;
            to = line.find_last_of("\"");
            if(line.find("STORAGE_PATH=")==0){
                fStoragePath=line.substr(from,to-from);
            }
        }
        if(configFile.eof()){configFile.clear();}
        configFile.close();
    }
    else{cout<<"SERVER -- Unable to open config file"<<endl;}
    TThread::UnLock();
    
    //start communication on socket
    cout<<"Starting server's communication"<<endl;
    StartCommunication();
}

AliStorageServerThread::~AliStorageServerThread()
{
    cout<<"SERVER -- AliStorageServerThread destructor called";
    if (fDatabase) {delete fDatabase;}
    cout<<" --- OK"<<endl;
}

void AliStorageServerThread::StartCommunication()
{
    AliZMQManager *eventManager = AliZMQManager::GetInstance();
    storageSockets socket = SERVER_COMMUNICATION_REP;
    eventManager->CreateSocket(socket);
    
    struct serverRequestStruct *request;
    
    bool receiveStatus = false;
    bool sendStatus = false;
    
    while(1)
    {
        cout<<"Server waiting for requests"<<endl;
        do{ // try to receive requests until success
            receiveStatus = eventManager->Get(request,socket);
        } while(receiveStatus == false);
        
        cout<<"Server received request:"<<request->messageType<<endl;
        switch(request->messageType)
        {
            case REQUEST_LIST_EVENTS:
            {
                cout<<"SERVER -- received request for list of events"<<endl;
                vector<serverListStruct> result = fDatabase->GetList(request->list);
                cout<<"SERVER -- got list from database"<<endl;
                sendStatus = eventManager->Send(result,socket);
                if(sendStatus){cout<<"SERVER -- list was sent"<<endl;}
                else{cout<<"SERVER -- couldn't send list"<<endl;}
                break;
            }
            case REQUEST_GET_EVENT:
            {
                cout<<"get event"<<endl;
                AliESDEvent *event = fDatabase->GetEvent(request->event);
                sendStatus = eventManager->Send(event,socket);
                delete event;
                break;
            }
            case REQUEST_GET_NEXT_EVENT:
            {
                cout<<"NEXT EVENT request received"<<endl;
                AliESDEvent *event = fDatabase->GetNextEvent(request->event);
                sendStatus = eventManager->Send(event,socket);
                delete event;
                break;
            }
            case REQUEST_GET_PREV_EVENT:
            {
                cout<<"PREV request"<<endl;
                AliESDEvent *event = fDatabase->GetPrevEvent(request->event);
                sendStatus = eventManager->Send(event,socket);
                delete event;
                break;
            }
            case REQUEST_GET_LAST_EVENT:
            {
                cout<<"LAST request"<<endl;
                AliESDEvent *event = fDatabase->GetLastEvent();
                sendStatus = eventManager->Send(event,socket);
                delete event;
                break;
            }
            case REQUEST_GET_FIRST_EVENT:
            {
                cout<<"FIRST request"<<endl;
                AliESDEvent *event = fDatabase->GetFirstEvent();
                sendStatus = eventManager->Send(event,socket);
                delete event;
                break;
            }
            case REQUEST_MARK_EVENT:
            {
                cout<<"MARK request"<<endl;
                struct eventStruct *markData  = &(request->event);
                sendStatus = eventManager->Send(MarkEvent(*markData),socket);
                break;
            }
            default:
            {
                cout<<"SERVER -- unknown request message"<<endl;
                sendStatus = false;
                break;
            }
        }
        if(sendStatus == false)
        {
            eventManager->RecreateSocket(socket);// if couldn't send, recreate socket to be able to receive messages (currently socket is in SEND state)
        }
        delete request;
    }
}

bool AliStorageServerThread::MarkEvent(struct eventStruct event)
{
    string pathToFile = fDatabase->GetFilePath(event);
    TFile *tmpFile = new TFile(pathToFile.c_str(),"read");
    if(!tmpFile)
    {
        cout<<"SERVER -- couldn't open temp file"<<endl;
        return false;
    }
    
    tmpFile->cd(Form("run%d",event.runNumber));
    
    AliESDEvent *eventToMark = (AliESDEvent*)gDirectory->Get(Form("event%d",event.eventNumber));
    if(!eventToMark)
    {
        cout<<"SERVER -- couldn't find such event"<<endl;
        if(tmpFile){tmpFile->Close();delete tmpFile;}
        return false;
    }
    cout<<"SERVER -- Marking event:"<<eventToMark->GetEventNumberInFile()<<endl;
    
    TFile *permFile = new TFile(Form("%s/permEvents.root",fStoragePath.c_str()),"update");//open/create perm file
    
    if(!permFile)
    {
        cout<<"SERVER -- Couldn't open perm file"<<endl;
        if(tmpFile){tmpFile->Close();delete tmpFile;}
        if(eventToMark){delete eventToMark;}
        return false;
    }
    
    //create new directory for this run
    TDirectory *currentRun;
    if((currentRun = permFile->mkdir(Form("run%d",event.runNumber))))
    {
        cout<<"SERVER -- creating new directory for this run"<<endl;
        currentRun->cd();
    }
    else
    {
        cout<<"SERVER -- opening existing directory for this run"<<endl;
        permFile->cd(Form("run%d",event.runNumber));
    }
    
    //try to add record to the database
    if(!fDatabase->MarkEvent(event))
    {
        cout<<"SERVER -- could not mark event in the database"<<endl;
        if(tmpFile){delete tmpFile;}
        if(eventToMark){delete eventToMark;}
        if(permFile){delete permFile;}
        return false;
    }
    
    eventToMark->Write(Form("event%d",event.eventNumber));
    permFile->Close();
    tmpFile->Close();
    
    if(tmpFile){delete tmpFile;}
    if(eventToMark){delete eventToMark;}
    if(permFile){delete permFile;}
    //	if(currentRun)delete currentRun;//this line crashes if there is no permanent file yet
    return true;
}

