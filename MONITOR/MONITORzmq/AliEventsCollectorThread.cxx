#include "AliEventsCollectorThread.h"
#include "AliZMQManager.h"

#include <TSystemDirectory.h>

#include <iostream>
#include <fstream>

using namespace std;

AliEventsCollectorThread::AliEventsCollectorThread(AliStorageClientThread *onlineReconstructionManager) :
fManager(onlineReconstructionManager),
fCollectorThread(0),
fCurrentFile(0),
fDatabase(0)
{
  fDatabase = new AliStorageDatabase();
    
    CheckCurrentStorageSize();
    
    // start collecting events in a thread
    fCollectorThread = new TThread("fCollectorThread",Dispatch,(void*)this);
    fCollectorThread->Run();
}


AliEventsCollectorThread::~AliEventsCollectorThread()
{
    if(fCollectorThread){delete fCollectorThread;}
    
    if(fCurrentFile){
        fCurrentFile->Close();
        delete fCurrentFile;
    }
    if(fDatabase){delete fDatabase;}
    if(fManager){delete fManager;}
}

void AliEventsCollectorThread::Kill()
{
    if(fCollectorThread){
        fCollectorThread->Join();
        fCollectorThread->Kill();
    }
}

void AliEventsCollectorThread::CollectorHandle()
{
    AliZMQManager *eventManager = AliZMQManager::GetInstance();
    eventManager->CreateSocket(EVENTS_SERVER_SUB);
    
    int chunkNumber=0;
    int previousChunkNumber=-1;
    int eventsInChunk=0;
    int previousRunNumber=-1;
    AliESDEvent *event = NULL;
    vector<struct eventStruct> eventsToUpdate;
    struct eventStruct currentEvent;
    
    int receiveStatus = false;
    
    while(1)
    {
        cout<<"CLIENT -- waiting for event..."<<endl;
        receiveStatus = eventManager->Get(event,EVENTS_SERVER_SUB);
        
        if (receiveStatus == 0){ // timeout
            continue;
        }
        else if(receiveStatus == -1){ // error, socket closed
            break;
        }
        else if(event && receiveStatus)
        {
            cout<<"CLIENT -- received event"<<endl;
            fManager->fReceivingStatus=STATUS_OK;
            
            if(event->GetRunNumber() != previousRunNumber)//when new run starts
            {
                cout<<"CLIENT -- new run started"<<endl;
                previousRunNumber = event->GetRunNumber();
                gSystem->Exec(Form("mkdir -p %s/run%d",fManager->fStoragePath.c_str(),event->GetRunNumber()));
                chunkNumber=0;
                eventsInChunk=0;
                
                TSystemDirectory dir(Form("%s/run%d",fManager->fStoragePath.c_str(),event->GetRunNumber()),
                                     Form("%s/run%d",fManager->fStoragePath.c_str(),event->GetRunNumber()));
                TList *files = dir.GetListOfFiles();
                if (files)
                {
                    TSystemFile *file;
                    string fname;
                    TIter next(files);
                    
                    while ((file=(TSystemFile*)next()))
                    {
                        fname = file->GetName();
                        
                        if (!file->IsDirectory())
                        {
                            int from = fname.find("chunk")+5;
                            int to = fname.find(".root");
                            
                            int maxChunkNumber = atoi(fname.substr(from,to-from).c_str());
                            
                            if(maxChunkNumber > chunkNumber)
                            {
                                chunkNumber = maxChunkNumber;
                            }
                        }
                    }
                    chunkNumber++;
                }
            }
            
            cout<<"CLIENT -- Received data. Event:"<<event->GetEventNumberInFile()<<"\trun:"<<event->GetRunNumber()<<endl;
            
            if(chunkNumber != previousChunkNumber)//when new chunk needs to be created
            {
                if(fCurrentFile)
                {
                    fCurrentFile->Close();
                    delete fCurrentFile;
                    fCurrentFile=0;
                }
                for(unsigned int i=0;i<eventsToUpdate.size();i++)
                {
		  TThread::Lock();
                    fDatabase->UpdateEventPath(eventsToUpdate[i],Form("%s/run%d/chunk%d.root",
                                                    fManager->fStoragePath.c_str(),
                                                    event->GetRunNumber(),
                                                    chunkNumber-1));
		    TThread::UnLock();
                }
                eventsToUpdate.clear();
                
                CheckCurrentStorageSize();
                
                fCurrentFile = new TFile(Form("%s/run%d/chunk%d.root", fManager->fStoragePath.c_str(),event->GetRunNumber(),chunkNumber),"recreate");
                
                previousChunkNumber = chunkNumber;
            }
            
            //create new directory for this run
            TDirectory *currentRun;
            if((currentRun = fCurrentFile->mkdir(Form("run%d",event->GetRunNumber()))))
            {
                cout<<"CLIENT -- creating new directory for this run"<<endl;
                currentRun->cd();
            }
            else
            {
                cout<<"CLIENT -- opening existing directory for this run"<<endl;
                fCurrentFile->cd(Form("run%d",event->GetRunNumber()));
            }
            
            if(0 != event->Write(Form("event%d",event->GetEventNumberInFile())))
                //fCurrentFile->WriteObject(event,Form("event%d",event->GetEventNumberInFile())))//if event was written to file
            {
                eventsInChunk++;
                
                if(eventsInChunk == fManager->fNumberOfEventsInFile)//if max events number in file was reached
                {
                    chunkNumber++;
                    eventsInChunk=0;
                }
                
                if(fManager->fSavingStatus!=STATUS_OK){fManager->fSavingStatus=STATUS_OK;}
            }
            else if(fManager->fSavingStatus!=STATUS_ERROR){fManager->fSavingStatus=STATUS_ERROR;}
            
            // save to event file as well:
            TFile *eventFile = new TFile(Form("%s/run%d/event%d.root", fManager->fStoragePath.c_str(),event->GetRunNumber(),eventsInChunk),"recreate");
            
            if((currentRun = eventFile->mkdir(Form("run%d",event->GetRunNumber()))))
            {
                cout<<"CLIENT -- creating new directory for this run"<<endl;
                currentRun->cd();
            }
            else
            {
                cout<<"CLIENT -- opening existing directory for this run"<<endl;
                eventFile->cd(Form("run%d",event->GetRunNumber()));
            }
            
            if(0 == event->Write(Form("event%d",event->GetEventNumberInFile())) &&
               fManager->fSavingStatus!=STATUS_ERROR)
            {
                fManager->fSavingStatus=STATUS_ERROR;
            }
            else
            {
                eventFile->Close();
                delete eventFile;
		TThread::Lock();
                
                cout<<"COLLECTOR -- events mask:"<<event->GetTriggerMask()<<endl;
                cout<<"COLLECTOR -- events mask next 50:"<<event->GetTriggerMaskNext50()<<endl;
                
                fDatabase->InsertEvent(event->GetRunNumber(),
                                       event->GetEventNumberInFile(),
                                       (char*)event->GetBeamType(),
                                       event->GetMultiplicity()->GetNumberOfTracklets(),
                                       Form("%s/run%d/event%d.root",fManager->fStoragePath.c_str(),
                                            event->GetRunNumber(),
                                            eventsInChunk),
                                       event->GetTriggerMask(),
                                       event->GetTriggerMaskNext50()
                                       );
		TThread::UnLock();
                currentEvent.runNumber = event->GetRunNumber();
                currentEvent.eventNumber = event->GetEventNumberInFile();
                eventsToUpdate.push_back(currentEvent);
            }
            delete event;event=0;
        }
        else
        {
            cout<<"CLIENT -- ERROR -- NO DATA!"<<endl;
            if(fManager->fReceivingStatus!=STATUS_ERROR){fManager->fReceivingStatus=STATUS_ERROR;}
        }
    }
    if(event){delete event;}
}


Long64_t AliEventsCollectorThread::GetSizeOfAllChunks()
{
    Long64_t totalStorageSize = 0;
    
    TSystemDirectory dir(fManager->fStoragePath.c_str(),fManager->fStoragePath.c_str());
    TList *listOfDirectories = dir.GetListOfFiles();
    
    if (!listOfDirectories){
        cout<<"CLIENT -- Storage directory is empty"<<endl;
        return 0;
    }
    TIter nextDirectory(listOfDirectories);
    TSystemFile *runDirectory;
    string directoryName;
    
    while ((runDirectory=(TSystemFile*)nextDirectory()))
    {
        directoryName=runDirectory->GetName();
        if (runDirectory->IsDirectory() && directoryName.find("run")==0)
        {
            TSystemDirectory dirChunks(Form("%s/%s",fManager->fStoragePath.c_str(),directoryName.c_str()),Form("%s/%s",fManager->fStoragePath.c_str(),directoryName.c_str()));
            TList *listOfChunks = dirChunks.GetListOfFiles();
            
            if(listOfChunks)
            {
                TIter nextChunk(listOfChunks);
                TSystemFile *chunk;
                string chunkFileName;
                
                while((chunk=(TSystemFile*)nextChunk()))
                {
                    chunkFileName = chunk->GetName();
                    if(!chunk->IsDirectory() && chunkFileName.find("chunk")==0)
                    {
                        TFile *tmpFile = new TFile(Form("%s/%s/%s",fManager->fStoragePath.c_str(),directoryName.c_str(),chunkFileName.c_str()),"read");
                        if(tmpFile)
                        {
                            totalStorageSize+=tmpFile->GetSize();
                            tmpFile->Close();
                            delete tmpFile;
                        }
                    }
                }
                if(chunk){delete chunk;}
            }
            if(listOfChunks){delete listOfChunks;}
        }
    }

    if(listOfDirectories){delete listOfDirectories;}
    if(runDirectory){delete runDirectory;}
    
    printf("CLIENT -- Total storage size:%lld\t(%.2f MB)\n",totalStorageSize,(float)totalStorageSize/(1000.*1000.));
    
    return totalStorageSize;
}

void AliEventsCollectorThread::CheckCurrentStorageSize()
{
    fManager->fCurrentStorageSize=GetSizeOfAllChunks();
    
    if(fManager->fCurrentStorageSize >  (float)fManager->fStorageOccupationLevel/100. * fManager->fMaximumStorageSize)
    {
        while(GetSizeOfAllChunks() > (float)fManager->fRemoveEventsPercentage/100. * fManager->fMaximumStorageSize)
        {
	  TThread::Lock();
            struct eventStruct oldestEvent = fDatabase->GetOldestEvent();
            string oldestEventPath = fDatabase->GetFilePath(oldestEvent);
	    TThread::UnLock();
            //remove oldest event
            cout<<"CLIENT -- Removing old events:"<<oldestEventPath<<endl;
            gSystem->Exec(Form("rm -f %s",oldestEventPath.c_str()));
	    TThread::Lock();
            fDatabase->RemoveEventsWithPath(oldestEventPath);
	    TThread::UnLock();
        }
    }
}



