#include "AliStorageClientThread.h"
#include "AliMultiplicity.h"
#include "AliStorageTypes.h"
#include "AliSocket.h"
#include "AliStorageEventManager.h"

//#include "zhelpers.hpp"
#include "zmq.hpp"

#include <sstream>
#include <signal.h>
#include <fstream>
#include <iostream>

#include <TSystemDirectory.h>
#include <TThread.h>
#include <TFile.h>

using namespace std;
using namespace zmq;

bool gClientQuit = false;    // signal flag
void GotSignalClient(int){gClientQuit = true;}

AliStorageClientThread::AliStorageClientThread() :
fConnectionStatus(STATUS_WAITING),
fReceivingStatus(STATUS_WAITING),
fSavingStatus(STATUS_WAITING),
fEventServer(""),
fCommunicationThread(0),
fCurrentFile(0),
fDatabase(0),
fCurrentStorageSize(0),
fMaximumStorageSize(0),
fStoragePath(""),
fNumberOfEventsInFile(0),
fStorageOccupationLevel(0),
fRemoveEventsPercentage(0)
{	
	// make sure that when program is closed destructor will be called
	struct sigaction sa;
	memset(&sa,0,sizeof(sa));
	sa.sa_handler = GotSignalClient;
	sigfillset(&sa.sa_mask);
	sigaction(SIGINT,&sa,NULL);
	
	//load storage parameters from file
	TThread::Lock();
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
			else if(line.find("MAX_SIZE=")==0)
			{
				fMaximumStorageSize=atoi(line.substr(from,to-from).c_str());
			}
			else if(line.find("MAX_OCCUPATION=")==0)
			{
				fStorageOccupationLevel=atoi(line.substr(from,to-from).c_str());
			}
			else if(line.find("REMOVE_PERCENT=")==0)
			{
				fRemoveEventsPercentage=atoi(line.substr(from,to-from).c_str());
			}
			else if(line.find("EVENTS_IN_FILE=")==0)
			{
				fNumberOfEventsInFile=atoi(line.substr(from,to-from).c_str());
			}
            else if(line.find("EVENT_SERVER=")==0)
			{
				fEventServer=line.substr(from,to-from);
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
		cout<<"CLIENT -- Unable to open config file"<<endl;
	}
	//create directory for storage if it doesn't exist
	gSystem->Exec(Form("mkdir -p %s",fStoragePath.c_str()));

	//create database class
	fDatabase = new AliStorageDatabase();
	TThread::UnLock();

	//check current storage size
	fCurrentStorageSize = GetSizeOfAllChunks();

	//create two-way commynication thread
	fCommunicationThread = new TThread("fCommunicationThread",
					   AliStorageClientThread::CommunicationHandler,(void*)this);
	fCommunicationThread->Run();
}

AliStorageClientThread::~AliStorageClientThread()
{
	cout<<"CLIENT -- AliStorageClientThread destructor called";
	if(fCurrentFile)
	{
		fCurrentFile->Close();
		delete fCurrentFile;
	}
	if(fCommunicationThread){delete fCommunicationThread;}
	cout<<" --- OK"<<endl;
}

void* AliStorageClientThread::CommunicationHandler(void *arg)
{
	AliStorageClientThread *clientInstance = static_cast<AliStorageClientThread*>(arg);

//create socket for two-way communication
	context_t *context = new context_t();
        socket_t *socket = new socket_t(*context,ZMQ_REP);
	socket->bind(Form("tcp://*:%d",gClientCommunicationPort));	

	AliStorageEventManager *eventManager = new AliStorageEventManager();
	
	message_t *requestMessage = new message_t();
	struct clientRequestStruct *request = new struct clientRequestStruct;
	struct clientRequestStruct *response = new struct clientRequestStruct;

	cout<<"CLIENT -- Communication stated"<<endl;
	
	while(!gClientQuit)
	{
		socket->recv(requestMessage);
		request = static_cast<struct clientRequestStruct*>(requestMessage->data());
		switch(request->messageType)
		{
		case REQUEST_CONNECTION:
			eventManager->Send((long)clientInstance->fConnectionStatus,socket);
			break;
		case REQUEST_RECEIVING:
			eventManager->Send((long)clientInstance->fReceivingStatus,socket);
			break;
		case REQUEST_SAVING:
			eventManager->Send((long)clientInstance->fSavingStatus,socket);
			break;
		case REQUEST_CURRENT_SIZE:
			eventManager->Send((long)clientInstance->fCurrentStorageSize,socket);
			break;
		case REQUEST_GET_PARAMS:
			response->maxStorageSize = clientInstance->fMaximumStorageSize;
			response->maxOccupation = clientInstance->fStorageOccupationLevel;
			response->removeEvents = clientInstance->fRemoveEventsPercentage;
			response->eventsInChunk = clientInstance->fNumberOfEventsInFile;

			eventManager->Send(response,socket);
			break;
		case REQUEST_SET_PARAMS:
			clientInstance->SetStorageParams(request->maxStorageSize,
							 request->maxOccupation,
							 request->removeEvents,
							 request->eventsInChunk);

			clientInstance->fMaximumStorageSize = request->maxStorageSize;
			clientInstance->fStorageOccupationLevel = request->maxOccupation;
			clientInstance->fRemoveEventsPercentage = request->removeEvents;
			clientInstance->fNumberOfEventsInFile = request->eventsInChunk;

			eventManager->Send(true,socket);
			break;
		default:break;
		}
	}
	if(context){delete context;}
	if(socket){delete socket;}
	if(eventManager){delete eventManager;}
	return nullptr;
}

void AliStorageClientThread::SetStorageParams(int maxStorageSize,int maxOccupation,int removeEvents,int eventsInChunk)
{
	cout<<maxStorageSize<<endl<<maxOccupation<<endl<<removeEvents<<endl<<eventsInChunk<<endl;


	TThread::Lock();
	ifstream configFile (Form("%s/STORAGE/setupStorageDatabase.sh",
				  gSystem->Getenv("ALICE_ROOT")));
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
			if(line.find("MAX_SIZE=")==0)
			{
				tmpLine = Form("MAX_SIZE=\"%d\"",maxStorageSize);
			}
			else if(line.find("MAX_OCCUPATION=")==0)
			{
				tmpLine = Form("MAX_OCCUPATION=\"%d\"",maxOccupation);
			}
			else if(line.find("REMOVE_PERCENT=")==0)
			{
				tmpLine = Form("REMOVE_PERCENT=\"%d\"",removeEvents);
			}
			else if(line.find("EVENTS_IN_FILE=")==0)
			{
				tmpLine = Form("EVENTS_IN_FILE=\"%d\"",eventsInChunk);
			}
			tmpLine += "\n";
			tmpFile << tmpLine;
		}
		if(configFile.eof())
		{
			configFile.clear();
		}
		configFile.close();
		tmpFile.close();
		rename("tmpFile.bla",Form("%s/STORAGE/setupStorageDatabase.sh",
					  gSystem->Getenv("ALICE_ROOT")));
	}
	else
	{
		cout<<"CLIENT -- Unable to open config file"<<endl;
	}
	TThread::UnLock();
}

Long64_t AliStorageClientThread::GetSizeOfAllChunks()
{
	Long64_t totalStorageSize = 0;

	TSystemDirectory dir(fStoragePath.c_str(),fStoragePath.c_str());
	TList *listOfDirectories = dir.GetListOfFiles();

	if (!listOfDirectories)
	{
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
			TSystemDirectory dirChunks(Form("%s/%s",fStoragePath.c_str(),directoryName.c_str()),Form("%s/%s",fStoragePath.c_str(),directoryName.c_str()));
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
						TFile *tmpFile = new TFile(Form("%s/%s/%s",fStoragePath.c_str(),directoryName.c_str(),chunkFileName.c_str()),"read");

						totalStorageSize+=tmpFile->GetSize();
						tmpFile->Close();
						if(tmpFile){delete tmpFile;}
					}
				}
				if(chunk){delete chunk;}
			}
			if(listOfChunks){delete listOfChunks;}
		}
	}

	//tmpFiles.clear();
	if(listOfDirectories){delete listOfDirectories;}
	if(runDirectory){delete runDirectory;}
	
	printf("CLIENT -- Total storage size:%lld\t(%.2f MB)\n",totalStorageSize,(float)totalStorageSize/(1000.*1000.));

	return totalStorageSize;
}

void AliStorageClientThread::CollectData()
{
	//connect with events' sockets
	context_t *context = new context_t(1);
	socket_t *socket = new socket_t(*context,ZMQ_SUB);
	socket->setsockopt(ZMQ_SUBSCRIBE,"",0);
	socket->connect(Form("tcp://%s:%d",fEventServer.c_str(),gEventsSubscriberPort));

	if(socket)
	{
		cout<<"CLIENT -- Successfully connected to events' server"<<endl;
		fConnectionStatus = STATUS_OK;
	}
       	else
	{
		cout<<"CLIENT -- ERROR - could not connect to events' server"<<endl;
		fConnectionStatus = STATUS_ERROR;
	}

	AliStorageEventManager *eventManager = new AliStorageEventManager();
	
	int chunkNumber=0;
	int previousChunkNumber=-1;
	int eventsInChunk=0;
	int previousRunNumber=-1;
	AliESDEvent *event = NULL;

	while(!gClientQuit)
	{		
		event = eventManager->GetEvent(socket);

		if(event)
		{
	       		fReceivingStatus=STATUS_OK;
			
			if(event->GetRunNumber() != previousRunNumber)//when new run starts
			{
				cout<<"CLIENT -- new run started"<<endl;
				previousRunNumber = event->GetRunNumber();
				gSystem->Exec(Form("mkdir -p %s/run%d",fStoragePath.c_str(),event->GetRunNumber()));
				chunkNumber=0;
				eventsInChunk=0;
				
				TSystemDirectory dir(Form("%s/run%d",fStoragePath.c_str(),event->GetRunNumber()),
						     Form("%s/run%d",fStoragePath.c_str(),event->GetRunNumber()));
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
				fCurrentStorageSize=GetSizeOfAllChunks();
				CheckCurrentStorageSize();
				
				fCurrentFile = new TFile(Form("%s/run%d/chunk%d.root", fStoragePath.c_str(),event->GetRunNumber(),chunkNumber),"recreate");

				previousChunkNumber = chunkNumber;
			}
			
			if(0 != fCurrentFile->WriteObject(event,Form("event%d",event->GetEventNumberInFile())))//if event was written to file
			{
				fDatabase->InsertEvent(event->GetRunNumber(),
					      event->GetEventNumberInFile(),
					      (char*)event->GetBeamType(),
						       event->GetMultiplicity()->GetNumberOfTracklets(),Form("%s/run%d/chunk%d.root",fStoragePath.c_str(),event->GetRunNumber(),chunkNumber));
				
				eventsInChunk++;
				
				if(eventsInChunk == fNumberOfEventsInFile)//if max events number in file was reached
				{
					chunkNumber++;
					eventsInChunk=0;
				}
				
				if(fSavingStatus!=STATUS_OK)
				{
					fSavingStatus=STATUS_OK;
				}
			}
			else if(fSavingStatus!=STATUS_ERROR)
			{
				fSavingStatus=STATUS_ERROR;
			}
			delete event;event=0;
		}
		else if(fReceivingStatus!=STATUS_ERROR)
		{
			cout<<"CLIENT -- ERROR -- NO DATA!"<<endl;
			fReceivingStatus=STATUS_ERROR;
		}
	}
	if(event){delete event;}
}


void AliStorageClientThread::CheckCurrentStorageSize()
{
	if(fCurrentStorageSize >  (float)fStorageOccupationLevel/100. * fMaximumStorageSize)
	{
		while(GetSizeOfAllChunks() > (float)fRemoveEventsPercentage/100. * fMaximumStorageSize)
		{
			struct eventStruct oldestEvent = fDatabase->GetOldestEvent();
			string oldestEventPath = fDatabase->GetFilePath(oldestEvent);
			//remove oldest event
			cout<<"CLIENT -- Removing old events:"<<oldestEventPath<<endl;
			gSystem->Exec(Form("rm -f %s",oldestEventPath.c_str()));
			fDatabase->RemoveEvent(fDatabase->GetOldestEvent());
		}
	}
}
