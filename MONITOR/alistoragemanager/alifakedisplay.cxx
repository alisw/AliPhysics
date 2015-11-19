#include "AliZMQManager.h"
#include "AliStorageTypes.h"
#include "AliESDEvent.h"
#include "AliESDRun.h"

#include <iostream>

#include "zmq.h"

#include <TThread.h>
#include <TFile.h>
#include "TMutex.h"
#include <TTimeStamp.h>

using namespace std;

// global variables:
AliESDEvent *currentEvent[2];
TTree *currentTree[2];
TMutex myMutex;
int eventInUse=1;
int writingToEventIndex=0;
bool isNewEventAvaliable=false;
//-----------------

static void* GetNextEvent(void*);

int main(int argc, char **argv)
{
    if(argc<2)
    {
        cout<<"Usage: alifakedisplay <mode>"<<endl;
        cout<<"mode:"<<endl;
        cout<<"0 - connect directly to reconstruction socket"<<endl;
        cout<<"1 - connect to Storage Manager and ask for last event"<<endl;
        cout<<"2 - connect to alieventserver, receive in thread and use in main thread"<<endl;
        cout<<"3 - connect to Storage Manager, download list of perm events and ask for them in loop"<<endl;
        return 0;
    }
    
    storageSockets socket;
    AliZMQManager *manager = AliZMQManager::GetInstance();
    AliESDEvent *event;
    struct recPointsStruct *files;
    
    if(atoi(argv[1])==0)
    {
        manager->CreateSocket(EVENTS_SERVER_SUB);
        
        while(1)
        {
            cout<<"waiting for event..."<<flush;
            manager->Get(event,EVENTS_SERVER_SUB);
            
            if(event)
            {
                if(event->GetRunNumber() > 0)
                {
                    cout<<"Received event. Run:"<<event->GetRunNumber()<<"\t event:"<<event->GetEventNumberInFile()<<endl;
                    
                    //                TTimeStamp ts(event->GetTimeStamp());
                    //                cout<<"Timestamp:"<<ts.AsString()<<endl;
                    
                    //
                    //                cout<<event->GetPeriodNumber()<<endl;
                    //                cout<<event->GetOrbitNumber()<<endl;
                    //                cout<<event->GetBunchCrossNumber()<<endl;
                    for(int i=0;i<100;i++)
                    {
                        if(strcmp(event->GetESDRun()->GetTriggerClass(i),""))
                        {
                            cout<<event->GetESDRun()->GetTriggerClass(i)<<endl;
                        }
                    }
                }
                delete event;event=0;
            }
            else
            {
                cout<<"NO EVENT"<<endl;
            }
            
            // Receiving RecPoints:
            /*
             
             files = manager->GetFiles(ITS_POINTS_SUB);
             
             
             cout<<"Received file."<<endl;
             
             //file->Write("ITS.RecPoints.root",kOverwrite);
             if(files->files[0]){files->files[0]->SaveAs("ITS.RecPoints.root","recreate");}
             if(files->files[1]){files->files[1]->SaveAs("TOF.RecPoints.root","recreate");}
             if(files->files[2]){files->files[2]->SaveAs("galice.root","recreate");}
             
             for(int i=0;i<10;i++)
             {
             if(files->files[i])
             {
             files->files[i]->Close();
             delete files->files[i];files->files[i]=0;
             }
             }
             sleep(2);
             */
        }
    }
    else if(atoi(argv[1])==1)
    {
        socket = SERVER_COMMUNICATION_REQ;
        manager->CreateSocket(socket);
        
        struct serverRequestStruct *requestMessage = new struct serverRequestStruct;
        requestMessage->messageType = REQUEST_GET_LAST_EVENT;
        
        while(1)
        {
            sleep(1);
            if(!manager->Send(requestMessage,socket))
            {
                cout<<"Couldn't send request"<<endl;
                continue;
            }
            
            if(!manager->Get(event,socket))
            {
                cout<<"Couldn't receive response"<<endl;
                manager->RecreateSocket(socket);
                continue;
            }
            
            if(event)
            {
                cout<<"Last event - Run:"<<event->GetRunNumber()<<"\t event:"<<event->GetEventNumberInFile()<<endl;
                delete event;
            }
            else
            {
                cout<<"NO EVENT"<<endl;
            }
        }
    }
    else if(atoi(argv[1])==2)
    {
        TThread *getEventThread = new TThread("getEventThread",GetNextEvent,nullptr);
        getEventThread->Run();
        
        int counter=0;
        
        while(1)
        {
            counter++;
            if(isNewEventAvaliable)
            {
                cout<<"new event"<<endl;
                myMutex.Lock();
                if(writingToEventIndex == 0) eventInUse = 0;
                else if(writingToEventIndex == 1) eventInUse = 1;
                cout<<"Using:"<<eventInUse<<endl;
                
                if(currentEvent[eventInUse])
                {
                    if(currentEvent[eventInUse]->GetRunNumber() >= 0)
                    {
                        cout<<"CURRENT EVENT:"<<currentEvent[eventInUse]->GetEventNumberInFile()<<endl;
                    }
                }
                isNewEventAvaliable = false;
                myMutex.UnLock();
            }
            else{cout<<"No new event is avaliable."<<endl;}
            
            sleep(2);
        }
    }
    else if(atoi(argv[1])==3)
    {
        socket = SERVER_COMMUNICATION_REQ;
        
        struct serverRequestStruct *requestMessage = new struct serverRequestStruct;
        requestMessage->messageType = REQUEST_LIST_EVENTS;
        //        requestMessage->list = list;
        
        //        struct listRequestStruct list;
        requestMessage->runNumber[0]=0;
        requestMessage->runNumber[1]=999999;
        requestMessage->eventNumber[0]=0;
        requestMessage->eventNumber[1]=999999;
        requestMessage->marked[0]=1;
        requestMessage->marked[1]=1;
        requestMessage->multiplicity[0]=1;
        requestMessage->multiplicity[1]=999999;
        strcpy(requestMessage->system[0],"p-p");
        strcpy(requestMessage->system[1],"A-A");
        
        
        manager->Send(requestMessage,socket);
        vector<serverListStruct> *tmpVector;
        manager->Get(tmpVector,socket);
        vector<serverListStruct> &receivedList = *tmpVector;
        
        cout<<"FAKE DISPLAY -- received list of marked events"<<endl;
        
        for(int i=0;i<receivedList.size();i++){
            cout<<"ev:"<<receivedList[i].eventNumber<<endl;
        }
        
        int iter=0;
        while(1)
        {
            if(iter<receivedList.size())
            {
                //                struct eventStruct mark;
                requestMessage->eventsRunNumber = receivedList[iter].runNumber;
                requestMessage->eventsEventNumber = receivedList[iter].eventNumber;
                
                requestMessage->messageType = REQUEST_GET_EVENT;
                //                requestMessage->event = mark;
                
                manager->Send(requestMessage,socket);
                manager->Get(event,socket);
                
                if(event)
                {
                    cout<<"i:"<<iter<<"\tevent:"<<event->GetEventNumberInFile()<<endl;
                }
                else
                {
                    cout<<"Received NO event"<<endl;
                }
                iter++;
                sleep(1);
            }
            else{iter=0;}
        }
    }
    else if(atoi(argv[1])==4)
    {
        cout<<"***************************"<<endl;
        cout<<"simulation of server thread"<<endl;
        cout<<"***************************"<<endl;
        socket = SERVER_COMMUNICATION_REP;
        
        long request;
        manager->Get(&request,socket);
        cout<<"\nreceived request:"<<request<<endl;
        
        if(request==55)
        {
            vector<serverListStruct> eventsVector;
            
            for(int i=0;i<5;i++)
            {
                serverListStruct resultList;
                
                resultList.runNumber = i;
                resultList.eventNumber = i;
                strcpy(resultList.system, "pp");
                resultList.multiplicity = 77;
                resultList.marked = 0;
                
                eventsVector.push_back(resultList);
            }
            manager->Send(eventsVector,socket);
            cout<<"sent reply"<<endl;
        }
        
        cout<<"waiting for bool"<<endl;
        bool boolMessage;
        manager->Get(&boolMessage,socket);
        cout<<"received bool:"<<boolMessage<<endl;
        
        cout<<"sending server request struct"<<endl;
        struct serverRequestStruct *srs = new struct serverRequestStruct;
        
        srs->messageType = REQUEST_GET_EVENT;
        //        struct eventStruct es;
        srs->eventsRunNumber = 197669;
        srs->eventsEventNumber = 123;
        //        srs->event = es;
        
        manager->Send(srs,socket);
        cout<<"server request struct sent"<<endl;
        
        cout<<"waiting for client struct"<<endl;
        struct clientRequestStruct *crs;
        manager->Get(crs,socket);
        cout<<"received client request struct:"<<crs->messageType<<"\tmaxOcc:"<<crs->maxOccupation<<"\tmaxsize:"<<crs->maxStorageSize<<endl;
        
        
        cout<<"sending AliESDevent"<<endl;
        TFile *esdFile = TFile::Open("~/Desktop/AliESDs.root");
        AliESDEvent *event = new AliESDEvent();
        TTree *esdTree = (TTree*)esdFile->Get("esdTree");
        event->ReadFromTree(esdTree);
        
        esdTree->GetEvent(0);
        cout<<"run number:"<<event->GetRunNumber()<<endl;
        cout<<"tracks:"<<event->GetNumberOfTracks()<<endl;
        
        manager->Send(event,socket);
        cout<<"event sent"<<endl;
    }
    else if(atoi(argv[1])==5)
    {
        cout<<"****************"<<endl;
        cout<<"simulation of ED"<<endl;
        cout<<"****************"<<endl;
        socket = SERVER_COMMUNICATION_REQ;
        
        cout<<"\nsending long request"<<endl;
        manager->Send((long)55,socket);
        
        vector<serverListStruct> *tmpVector;
        manager->Get(tmpVector,socket);
        vector<serverListStruct> &result = *tmpVector;
        
        cout<<"received vector list:"<<endl;
        for(int i=0;i<result.size();i++)
        {
            cout<<result[i].runNumber<<"\t"<<result[i].system<<endl;
        }
        cout<<"sending bool message:"<<true<<endl;
        manager->Send(true,socket);
        cout<<"bool message sent"<<endl;
        
        cout<<"waiting for server request struct"<<endl;
        struct serverRequestStruct *srs;
        manager->Get(srs,socket);
        
        cout<<"received server request struct:"<<srs->messageType<<"\trun:"<<srs->eventsRunNumber<<"\tevent:"<<srs->eventsEventNumber<<endl;
        
        cout<<"sending client request struct"<<endl;
        struct clientRequestStruct *crs = new struct clientRequestStruct;
        
        crs->messageType = REQUEST_SAVING;
        crs->maxStorageSize = 234;
        crs->maxOccupation = 0;
        crs->removeEvents = 0;
        crs->eventsInChunk = 0;
        
        manager->Send(crs,socket);
        cout<<"client request struct sent"<<endl;
        
        cout<<"waiting for AliESDevent"<<endl;
        AliESDEvent *event = NULL;
        manager->Get(event,socket);
        
        cout<<"received event with run number"<<event->GetRunNumber()<<"\ttracks:"<<event->GetNumberOfTracks()<<endl;
    }
    
    return 0;
}

void* GetNextEvent(void*)
{
    AliZMQManager *eventManager = AliZMQManager::GetInstance();
    
    currentEvent[0]=0;
    currentEvent[1]=0;
    currentTree[0]=0;
    currentTree[1]=0;
    AliESDEvent *tmpEvent;
    
    while(1)
    {
        //if(tmpEvent){delete tmpEvent;tmpEvent=0;}
        eventManager->Get(tmpEvent,EVENTS_SERVER_SUB);
        
        if(tmpEvent)
        {
            if(tmpEvent->GetRunNumber()>=0)
            {
                myMutex.Lock();
                if(eventInUse == 0){writingToEventIndex = 1;}
                else if(eventInUse == 1){writingToEventIndex = 0;}
                
                if(currentEvent[writingToEventIndex])
                {
                    delete currentEvent[writingToEventIndex];
                    currentEvent[writingToEventIndex]=0;
                }
                currentEvent[writingToEventIndex] = tmpEvent;
                isNewEventAvaliable = true;
                myMutex.UnLock();
            }
        }	
    }
}
