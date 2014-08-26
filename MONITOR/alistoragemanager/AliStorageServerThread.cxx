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
	TThread::Lock();
	fDatabase = new AliStorageDatabase();
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
			if(line.find("STORAGE_PATH=")==0)
			{
				fStoragePath=line.substr(from,to-from);
			}
		}
		if(configFile.eof())
		{
			configFile.clear();
		}
		configFile.close();
	}
	else
	{
		cout<<"SERVER -- Unable to open config file"<<endl;
	}
	TThread::UnLock();

	//start communication on socket
	StartCommunication();
}

AliStorageServerThread::~AliStorageServerThread()
{
	cout<<"SERVER -- AliStorageServerThread destructor called";	
	cout<<" --- OK"<<endl;
}

void AliStorageServerThread::StartCommunication()
{
	AliStorageEventManager *eventManager = AliStorageEventManager::GetEventManagerInstance();
	storageSockets socket = SERVER_COMMUNICATION_REP;
	eventManager->CreateSocket(socket);

	struct serverRequestStruct *request;
	
	while(1)
	{
		request = eventManager->GetServerStruct(socket);
		
		switch(request->messageType)
		{
		case REQUEST_LIST_EVENTS:
		{
			vector<serverListStruct> result = fDatabase->GetList(request->list);
			eventManager->Send(result,socket);
			break;
		}
		case REQUEST_GET_EVENT:
		{
			AliESDEvent *event = fDatabase->GetEvent(request->event);
			eventManager->Send(event,socket);
			delete event;
		    	break;
		}
		case REQUEST_GET_NEXT_EVENT:
		{
			AliESDEvent *event = fDatabase->GetNextEvent(request->event);
			eventManager->Send(event,socket);
			delete event;
			break;
		}
		case REQUEST_GET_LAST_EVENT:
		{
			AliESDEvent *event = fDatabase->GetLastEvent();
			eventManager->Send(event,socket);
			delete event;
			break;
		}
		case REQUEST_MARK_EVENT:
		{
			struct eventStruct *markData  = &(request->event);
			eventManager->Send(MarkEvent(*markData),socket);
			break;
		}
		default:break;
		}
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
	AliESDEvent *eventToMark = (AliESDEvent*)tmpFile->Get(Form("event%d",event.eventNumber));
	if(!eventToMark)
	{
		cout<<"SERVER -- couldn't find such event"<<endl;
		if(tmpFile){delete tmpFile;}
		return false;
	}
	cout<<"SERVER -- Marking event:"<<eventToMark->GetEventNumberInFile()<<endl;
		
	TFile *permFile = new TFile(Form("%s/permEvents.root",fStoragePath.c_str()),"update");//open/create perm file
	
	if(!permFile)
	{
		cout<<"SERVER -- Couldn't open perm file"<<endl;
		if(tmpFile){delete tmpFile;}
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

