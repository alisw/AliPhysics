#include "AliCommunicationThread.h"
#include "AliStorageEventManager.h"

#include <iostream>
#include <fstream>

using namespace std;

AliCommunicationThread::AliCommunicationThread(AliStorageClientThread *onlineReconstructionManager) :
    fFinished(false),
    fManager(onlineReconstructionManager),
    fCommunicationThread(0)
{
    //create two-way communication thread
    fCommunicationThread = new TThread("fCommunicationThread",Dispatch,(void*)this);
    fCommunicationThread->Run();
}

AliCommunicationThread::~AliCommunicationThread()
{
    if(fCommunicationThread){delete fCommunicationThread;}
}

void AliCommunicationThread::Kill()
{
    if(fCommunicationThread)
    {
        fFinished=true;
        fCommunicationThread->Join();
        fCommunicationThread->Kill();
    }
}

void AliCommunicationThread::CommunicationHandle()
{
    AliStorageEventManager *eventManager = AliStorageEventManager::GetEventManagerInstance();
    storageSockets socket = CLIENT_COMMUNICATION_REP;
    eventManager->CreateSocket(socket);
    
    struct clientRequestStruct *request;
    struct clientRequestStruct *response = new struct clientRequestStruct;
    
    cout<<"CLIENT -- Communication stated"<<endl;
    
    // mutex mtx;
    
    while(!fFinished)
    {
        cout<<"COMMUNICATION -- waiting for requests"<<endl;
        request = eventManager->GetClientStruct(socket,5000);
        
        if(request)
        {
	  //      lock_guard<mutex> lock(mtx);
            cout<<"COMMUNICATION -- received request"<<endl;
            switch(request->messageType)
            {
                case REQUEST_CONNECTION:
                    eventManager->Send((long)fManager->fConnectionStatus,socket);
                    break;
                case REQUEST_RECEIVING:
                    eventManager->Send((long)fManager->fReceivingStatus,socket);
                    break;
                case REQUEST_SAVING:
                    eventManager->Send((long)fManager->fSavingStatus,socket);
                    break;
                case REQUEST_CURRENT_SIZE:
                    eventManager->Send((long)fManager->fCurrentStorageSize,socket);
                    break;
                case REQUEST_GET_PARAMS:
                    response->maxStorageSize = fManager->fMaximumStorageSize;
                    response->maxOccupation = fManager->fStorageOccupationLevel;
                    response->removeEvents = fManager->fRemoveEventsPercentage;
                    response->eventsInChunk = fManager->fNumberOfEventsInFile;
                    
                    eventManager->Send(response,socket);
                    break;
                case REQUEST_SET_PARAMS:
                    SetStorageParams(request->maxStorageSize,
                                     request->maxOccupation,
                                     request->removeEvents,
                                     request->eventsInChunk);
                    
                    fManager->fMaximumStorageSize = request->maxStorageSize;
                    fManager->fStorageOccupationLevel = request->maxOccupation;
                    fManager->fRemoveEventsPercentage = request->removeEvents;
                    fManager->fNumberOfEventsInFile = request->eventsInChunk;
                    
                    eventManager->Send(true,socket);
                    break;
                default:break;
            }
            delete request;
        }
        else{cout<<"COMMUNICATION -- received NO request"<<endl;}
    }
}

void AliCommunicationThread::SetStorageParams(int maxStorageSize,int maxOccupation,int removeEvents,int eventsInChunk)
{
    cout<<maxStorageSize<<endl<<maxOccupation<<endl<<removeEvents<<endl<<eventsInChunk<<endl;
    
    TThread::Lock();
    ifstream configFile (GetConfigFilePath());
    ofstream tmpFile("tmpFile.bla");
    
    if (configFile.is_open())
    {
        string line;
        string tmpLine;
        int from,to;
        while(configFile.good())
        {
            getline(configFile,line);
            from = line.find("\"")+1;
            to = line.find_last_of("\"");
            tmpLine = line;
            if(line.find("MAX_SIZE=")==0){
                tmpLine = Form("MAX_SIZE=\"%d\"",maxStorageSize);
            }
            else if(line.find("MAX_OCCUPATION=")==0){
                tmpLine = Form("MAX_OCCUPATION=\"%d\"",maxOccupation);
            }
            else if(line.find("REMOVE_PERCENT=")==0){
                tmpLine = Form("REMOVE_PERCENT=\"%d\"",removeEvents);
            }
            else if(line.find("EVENTS_IN_FILE=")==0){
                tmpLine = Form("EVENTS_IN_FILE=\"%d\"",eventsInChunk);
            }
            tmpLine += "\n";
            tmpFile << tmpLine;
        }
        if(configFile.eof()){configFile.clear();}
        configFile.close();
        tmpFile.close();
        rename("tmpFile.bla",GetConfigFilePath());
    }
    else{cout<<"CLIENT -- Unable to open config file"<<endl;}
    TThread::UnLock();
}
