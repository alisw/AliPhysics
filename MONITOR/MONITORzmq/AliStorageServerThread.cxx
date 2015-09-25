#include "AliStorageServerThread.h"
#include "AliStorageTypes.h"
#include "AliESDEvent.h"
#include "AliCDBManager.h"
#include "AliCDBEntry.h"
#include "AliCDBPath.h"
#include "AliOnlineReconstructionUtil.h"
#include "AliTriggerConfiguration.h"
#include "AliTriggerClass.h"

#include <iostream>
#include <fstream>

#include <TEnv.h>
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
    
    int receiveStatus = false;
    int sendStatus = false;
    
    while(1)
    {
        cout<<"Server waiting for requests"<<endl;

        receiveStatus = eventManager->Get(request,socket);
        
        if(receiveStatus == 0){ // timeout
            continue;
        }
        else if(receiveStatus == -1){ // error, socket closed
            break;
        }
        
        cout<<"Server received request:"<<request->messageType<<endl;
        switch(request->messageType)
        {
            case REQUEST_LIST_EVENTS:
            {
                cout<<"SERVER -- received request for list of events"<<endl;
                struct listRequestStruct list;
                list.runNumber[0] = request->runNumber[0];
                list.runNumber[1] = request->runNumber[1];
                list.eventNumber[0] = request->eventNumber[0];
                list.eventNumber[1] = request->eventNumber[1];
                list.marked[0] = request->marked[0];
                list.marked[1] = request->marked[1];
                list.multiplicity[0] = request->multiplicity[0];
                list.multiplicity[1] = request->multiplicity[1];
                strcpy(list.system[0],request->system[0]);
                strcpy(list.system[1],request->system[1]);
                strcpy(list.triggerClass,request->triggerClass);
                
                vector<serverListStruct> result = fDatabase->GetList(list);
                cout<<"SERVER -- got list from database"<<endl;
                sendStatus = eventManager->Send(result,socket);
                if(sendStatus){cout<<"SERVER -- list was sent"<<endl;}
                else{cout<<"SERVER -- couldn't send list"<<endl;}
                break;
            }
            case REQUEST_GET_EVENT:
            {
                cout<<"get event"<<endl;
                struct eventStruct ev;
                ev.runNumber = request->eventsRunNumber;
                ev.eventNumber = request->eventsEventNumber;
                
                AliESDEvent *event = fDatabase->GetEvent(ev);
                sendStatus = eventManager->Send(event,socket);
                delete event;
                break;
            }
            case REQUEST_GET_NEXT_EVENT:
            {
                cout<<"NEXT EVENT request received"<<endl;
                
                struct eventStruct ev;
                ev.runNumber = request->eventsRunNumber;
                ev.eventNumber = request->eventsEventNumber;
                
                AliESDEvent *event = fDatabase->GetNextEvent(ev);
                sendStatus = eventManager->Send(event,socket);
                delete event;
                break;
            }
            case REQUEST_GET_PREV_EVENT:
            {
                cout<<"PREV request"<<endl;
                
                struct eventStruct ev;
                ev.runNumber = request->eventsRunNumber;
                ev.eventNumber = request->eventsEventNumber;
                
                AliESDEvent *event = fDatabase->GetPrevEvent(ev);
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
                struct eventStruct ev;
                ev.runNumber = request->eventsRunNumber;
                ev.eventNumber = request->eventsEventNumber;
                
                struct eventStruct *markData  = &(ev);
                sendStatus = eventManager->Send(MarkEvent(*markData),socket);
                break;
            }
            case REQUEST_GET_TRIGGER_LIST:
            {
                cout<<"Get trigger list request"<<endl;
                vector<string100> classes = GetTriggerClasses();
                cout<<"SERVER -- got set from database"<<endl;
                sendStatus = eventManager->Send(classes,socket);
                break;
            }
            default:
            {
                cout<<"SERVER -- unknown request message"<<endl;
                sendStatus = false;
                break;
            }
        }
        if(sendStatus == 0) // timeout
        {
            eventManager->RecreateSocket(socket);// if couldn't send, recreate socket to be able to receive messages (currently socket is in SEND state)
        }
        else if(sendStatus == -1){ // error, socket closed
            break;
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

vector<string100> AliStorageServerThread::GetTriggerClasses()
{
    set<string> setResult;
    
    // craete CDB manager:
    TEnv settings;
    settings.ReadFile(AliOnlineReconstructionUtil::GetPathToServerConf(), kEnvUser);
    const char *cdbPath = settings.GetValue("cdb.defaultStorage", "");
    AliCDBManager *man = AliCDBManager::Instance();
    man->SetDefaultStorage(cdbPath);
    
    // get list of runs stored in Storage's database:
    vector<int> runs = fDatabase->GetListOfRuns();
    
    // get trigger classes for all runs:
    AliCDBEntry *cdbEntry;
    AliTriggerConfiguration *cfg;
    TObjArray trarr;
    AliCDBPath path("GRP/CTP/Config");

    AliTriggerClass* trgclass;
    
    for(int i=0;i<runs.size();i++)
    {
        man->SetRun(runs[i]);
        cdbEntry = man->Get(path);
        cfg = (AliTriggerConfiguration*)cdbEntry->GetObject();
        trarr = cfg->GetClasses();
        
        for (int j=0;j<trarr.GetEntriesFast();j++)
        {
            trgclass = (AliTriggerClass*)trarr.At(j);
            if(strcmp(trgclass->GetName(),"NONE")!=0)
            {
                setResult.insert(trgclass->GetName());
            }
            else
            {
                cout<<"Trigger classes not defined for run "<<runs[i]<<endl;
                cout<<"This may indicate problems with OCDB"<<endl;
            }
        }
    }
    
    vector<string100> result;
    set<string>::iterator it;
    
    for (it = setResult.begin(); it != setResult.end(); ++it)
    {
        string str = *it;
        string100 cls(str.c_str());
        result.push_back(cls);
    }
    return result;
}







