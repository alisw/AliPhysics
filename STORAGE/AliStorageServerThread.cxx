#include "AliStorageServerThread.h"
#include "AliStorageTypes.h"
#include "AliESDEvent.h"

#include "zmq.hpp"
#include <iostream>
#include <fstream>

#include <TFile.h>
#include <TThread.h>

using namespace std;
using namespace zmq;

AliStorageServerThread::AliStorageServerThread() :
	fEventManager(0),
	fDatabase(0),
	fStoragePath("")
{
	TThread::Lock();
	fDatabase = new AliStorageDatabase();
	//load parameters from config file
      	ifstream configFile (Form("%s/STORAGE/setupStorageDatabase.sh",
				  gSystem->Getenv("ALICE_ROOT")));

	
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
	fEventManager = new AliStorageEventManager();
	StartCommunication();
}

AliStorageServerThread::~AliStorageServerThread()
{
	cout<<"SERVER -- AliStorageServerThread destructor called";
	
	cout<<" --- OK"<<endl;
}

void AliStorageServerThread::StartCommunication()
{
	//create two-way communication socket
	context_t *context = new context_t(1);
	socket_t *socket = new socket_t(*context,ZMQ_REP);
	socket->bind(Form("tcp://*:%d",gServerCommunicationPort));	
	
	message_t *request = new message_t;
	message_t *reply;
	struct serverRequestStruct *requestMessage = new struct serverRequestStruct;
	char *buffer;
	
	while(1)
	{
		socket->recv(request);
		requestMessage = static_cast<struct serverRequestStruct*>(request->data());
		
		switch(requestMessage->messageType)
		{
		case REQUEST_LIST_EVENTS:
		{
			vector<serverListStruct> result = fDatabase->GetList(requestMessage->list);
			fEventManager->Send(result,socket);
			break;
		}
		case REQUEST_GET_EVENT:
		{
			AliESDEvent *event = fDatabase->GetEvent(requestMessage->event);
			fEventManager->Send(event,socket);
			delete event;
		    	break;
		}
		case REQUEST_GET_NEXT_EVENT:
		{
			AliESDEvent *event = fDatabase->GetNextEvent(requestMessage->event);
			fEventManager->Send(event,socket);
			delete event;
			break;
		}
		case REQUEST_GET_LAST_EVENT:
		{
			AliESDEvent *event = fDatabase->GetLastEvent();
			fEventManager->Send(event,socket);
			delete event;
			break;
		}
		case REQUEST_MARK_EVENT:
		{
			struct eventStruct *markData  = &(requestMessage->event);
			buffer =(char*)(MarkEvent(*markData) ? "true" : "false");
			reply = new message_t((void*)buffer,sizeof(char*),0);
			socket->send(*reply);
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

