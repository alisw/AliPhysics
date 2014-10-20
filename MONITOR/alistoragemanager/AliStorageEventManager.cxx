#include "AliStorageEventManager.h"

#include <iostream>
#include <sstream>
#include <fstream>

#include <TList.h>
#include <TStreamerInfo.h>
#include <TThread.h>

#include "zmq.hpp"

#include "AliESDEvent.h"
#include "AliESDtrack.h"
#include "AliTrackPointArray.h"
#include "AliESDfriendTrack.h"
#include "AliExternalTrackParam.h"
#include "AliTrackerBase.h"
#include "AliTracker.h"

using namespace zmq;
using namespace std;

AliStorageEventManager *AliStorageEventManager::fManagerInstance = 0;

AliStorageEventManager::AliStorageEventManager()
{
	//read config file
	TThread::Lock();
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
			if(line.find("STORAGE_SERVER=")==0)
			{
				fStorageServer=line.substr(from,to-from);
			}
			else if(line.find("EVENT_SERVER=")==0)
			{
				fEventServer=line.substr(from,to-from);
			}
			else if(line.find("STORAGE_SERVER_PORT=")==0)
			{
				fStorageServerPort=atoi(line.substr(from,to-from).c_str());
			}
			else if(line.find("EVENT_SERVER_PORT=")==0)
			{
				fEventServerPort=atoi(line.substr(from,to-from).c_str());
			}
			else if(line.find("STORAGE_CLIENT_PORT=")==0)
			{
				fStorageClientPort=atoi(line.substr(from,to-from).c_str());
			}
			else if(line.find("XML_SERVER_PORT=")==0)
			{
				fXmlServerPort=atoi(line.substr(from,to-from).c_str());
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
		cout<<"EVENT MANAGER -- Unable to open config file"<<endl;
	}
	TThread::UnLock();
	
	for(int i=0;i<NUMBER_OF_SOCKETS;i++)
	{
		fContexts[i] = new context_t();
	}
}
AliStorageEventManager::~AliStorageEventManager()
{
	if(fManagerInstance){delete fManagerInstance;fManagerInstance=0;}
}

AliStorageEventManager* AliStorageEventManager::GetEventManagerInstance()
{
	TThread::Lock();
	if(fManagerInstance==0)
	{
		fManagerInstance = new AliStorageEventManager();
	}
	TThread::UnLock();
	return fManagerInstance;
}


void freeBuff (void *data, void *hint)
{
 //   free(data);
}

bool AliStorageEventManager::CreateSocket(storageSockets socket)
{
    cout<<"Creating socket:"<<socket<<endl;
    
    switch(socket)
    {
        case SERVER_COMMUNICATION_REQ:
        {
            fSockets[SERVER_COMMUNICATION_REQ] =
            new socket_t(*fContexts[SERVER_COMMUNICATION_REQ],ZMQ_REQ);
            try
            {
                fSockets[SERVER_COMMUNICATION_REQ]->connect(Form("tcp://%s:%d",fStorageServer.c_str(),fStorageServerPort));
            }
            catch (const zmq::error_t& e)
            {
                cout<<"MANAGER -- "<<e.what()<<endl;
                return 0;
            }
        }
            break;
        case SERVER_COMMUNICATION_REP:
        {
            fSockets[SERVER_COMMUNICATION_REP] =
            new socket_t(*fContexts[SERVER_COMMUNICATION_REP],ZMQ_REP);
            try
            {
                fSockets[SERVER_COMMUNICATION_REP]->bind(Form("tcp://*:%d",fStorageServerPort));
            }
            catch (const zmq::error_t& e)
            {
                cout<<"MANAGER -- "<<e.what()<<endl;
                return 0;
            }
        }
            break;
        case CLIENT_COMMUNICATION_REQ:
        {
            fSockets[CLIENT_COMMUNICATION_REQ] =
            new socket_t(*fContexts[CLIENT_COMMUNICATION_REQ],ZMQ_REQ);
            try
            {
                fSockets[CLIENT_COMMUNICATION_REQ]->connect(Form("tcp://%s:%d",fStorageServer.c_str(), fStorageClientPort));
            }
            catch (const zmq::error_t& e)
            {
                cout<<"MANAGER -- "<<e.what()<<endl;
                return 0;
            }
        }
            break;
        case CLIENT_COMMUNICATION_REP:
        {
            fSockets[CLIENT_COMMUNICATION_REP] =
            new socket_t(*fContexts[CLIENT_COMMUNICATION_REP],ZMQ_REP);
            try
            {
                fSockets[CLIENT_COMMUNICATION_REP]->bind(Form("tcp://*:%d",fStorageClientPort));
            }
            catch (const zmq::error_t& e)
            {
                cout<<"MANAGER -- "<<e.what()<<endl;
                return 0;
            }
        }
            break;
        case EVENTS_SERVER_PUB:
        {
            fSockets[EVENTS_SERVER_PUB] =
            new socket_t(*fContexts[EVENTS_SERVER_PUB],ZMQ_PUB);
            try
            {
                fSockets[EVENTS_SERVER_PUB]->bind(Form("tcp://*:%d",fEventServerPort));
            }
            catch (const zmq::error_t& e)
            {
                cout<<"MANAGER -- "<<e.what()<<endl;
                return 0;
            }
        }
            break;
        case EVENTS_SERVER_SUB:
        {
            fSockets[EVENTS_SERVER_SUB] =
            new socket_t(*fContexts[EVENTS_SERVER_SUB],ZMQ_SUB);
            fSockets[EVENTS_SERVER_SUB]->setsockopt(ZMQ_SUBSCRIBE,"",0);
            try
            {
                fSockets[EVENTS_SERVER_SUB]->connect(Form("tcp://%s:%d",fEventServer.c_str(),fEventServerPort));
            }
            catch (const zmq::error_t& e)
            {
                cout<<"MANAGER -- "<<e.what()<<endl;
                return 0;
                
            }
        }
            break;
        case XML_PUB:
        {
            fSockets[XML_PUB] =
            new socket_t(*fContexts[XML_PUB],ZMQ_PUB);
            try
            {
                fSockets[XML_PUB]->bind(Form("tcp://*:%d",fXmlServerPort));
            }
            catch (const zmq::error_t& e)
            {
                cout<<"MANAGER -- "<<e.what()<<endl;
                return 0;
            }
        }
            break;
        default:break;
    }
    return 1;
}

void AliStorageEventManager::Send(vector<serverListStruct> list,storageSockets socket)
{
	//send size of the struct first
	int numberOfRecords = list.size();
	message_t message(20);
	snprintf ((char *)message.data(), 20 ,"%d",numberOfRecords);

	fSockets[socket]->send(message);
	if(numberOfRecords==0)return;
    message_t *tmpMessage = new message_t();
	fSockets[socket]->recv(tmpMessage);//empty message just to keep req-rep order
    
	// //prepare message with event's list
	// char *buffer = reinterpret_cast<char*> (&list[0]);
	// message_t *reply = new message_t((void*)buffer,
	// 		      sizeof(serverListStruct)*numberOfRecords,0);
	// fSockets[socket]->send(*reply);
	// if(reply){delete reply;}

	message_t reply(sizeof(serverListStruct)*numberOfRecords);
	memcpy(reply.data(), reinterpret_cast<const char*> (&list[0]), sizeof(serverListStruct)*numberOfRecords);

	fSockets[socket]->send(reply);
    if(tmpMessage){delete tmpMessage;}
}

void AliStorageEventManager::Send(struct serverRequestStruct *request,storageSockets socket)
{
	char *buffer = (char*)(request);
	message_t *requestMessage = new message_t((void*)buffer,
					   sizeof(struct serverRequestStruct)
					   +sizeof(struct listRequestStruct)
					   +sizeof(struct eventStruct),freeBuff);
	fSockets[socket]->send(*requestMessage);
}

bool AliStorageEventManager::Send(struct clientRequestStruct *request,storageSockets socket,int timeout)
{
	pollitem_t items[1] =  {{*fSockets[socket],0,ZMQ_POLLIN,0}} ;
	
	
	char *buffer = (char*)(request);
	message_t *requestMessage = new message_t((void*)buffer,
						  sizeof(struct clientRequestStruct),freeBuff);

	try
	{
		fSockets[socket]->send(*requestMessage);
	}
	catch (const zmq::error_t& e)
	{
		cout<<"MANAGER -- "<<e.what()<<endl;
		cout<<e.num()<<endl;
		if(fSockets[socket]){delete fSockets[socket];fSockets[socket]=0;}

		CreateSocket(socket);
		delete requestMessage;
		return 0;
		
	}	
	if(timeout>=0)
	{
		if(poll (&items[0], 1, timeout)==0)
		{
		  delete requestMessage;
			return 0;
		}
	}
	delete requestMessage;		
	return 1;
}

void AliStorageEventManager::Send(long message,storageSockets socket)
{
	stringstream streamBuffer;
	streamBuffer << message;
	string stringBuffer = streamBuffer.str();
	char *buffer = (char*)stringBuffer.c_str();
	message_t *replyMessage = new message_t((void*)buffer,sizeof(stringBuffer),freeBuff);
    
	fSockets[socket]->send(*replyMessage);
	delete replyMessage;
	streamBuffer.str(string());
	streamBuffer.clear();
}

void AliStorageEventManager::Send(bool message,storageSockets socket)
{
	char *buffer;
	if(message==true)
	{
		buffer = (char*)("true");
	}
	else
	{
		buffer = (char*)("false");
	}
	message_t *replyMessage = new message_t((void*)buffer,sizeof(buffer),freeBuff);
	fSockets[socket]->send(*replyMessage);
	delete replyMessage;
}

void AliStorageEventManager::Send(AliESDEvent *event, storageSockets socket)
{
	TMessage tmess(kMESS_OBJECT);
	tmess.Reset();
	tmess.WriteObject(event);
	TMessage::EnableSchemaEvolutionForAll(kTRUE);

	int bufsize = tmess.Length();
	char* buf = (char*) malloc(bufsize * sizeof(char));
	memcpy(buf, tmess.Buffer(), bufsize);
	
	message_t message((void*)buf, bufsize, freeBuff);
	fSockets[socket]->send(message);
}

void AliStorageEventManager::SendAsXml(AliESDEvent *event,storageSockets socket)
{
	cout<<"SENDING AS XML"<<endl;
	stringstream bufferStream;
	bufferStream << "<?xml version=\"1.0\" encoding=\"utf-8\"?>"<<endl;
	bufferStream << "<ESD xmlns:xsi=\"http://www.w3.org/2001/XMLSchema-instance\" xmlns:xsd=\"http://www.w3.org/2001/XMLSchema\">" << endl;
	
	for(int i=0;i<event->GetNumberOfTracks();i++)
	{
		AliESDtrack *track = event->GetTrack(i);
		bufferStream << "\t<track mass=\""<<track->GetMass()<<"\">" <<endl;
		const AliTrackPointArray *array = track->GetTrackPointArray();

		if(array)
		{
			const float *x = array->GetX();
			const float *y = array->GetY();
			const float *z = array->GetZ();
			int n = array->GetNPoints();

			for(int j=0;j<n;j++)
			{
				bufferStream <<"\t\t<point>"<<endl;
				bufferStream <<"\t\t\t<x>"<< x[j] <<"</x>"<<endl;
				bufferStream <<"\t\t\t<y>"<< y[j] <<"</y>"<<endl;
				bufferStream <<"\t\t\t<z>"<< z[j] <<"</z>"<<endl;
				bufferStream <<"\t\t</point>"<<endl;
			}
		}
		else cout<<"no array"<<endl;

		bufferStream << "\t</track>"<<endl;
	}

	bufferStream << "</ESD>"<<endl;

	string bufferString = bufferStream.str();
	message_t message(bufferString.size());
	memcpy (message.data(), bufferString.data(), bufferString.size());	
	
	fSockets[socket]->send(message);
	cout<<"xml sent"<<endl;
}

vector<serverListStruct> AliStorageEventManager::GetServerListVector(storageSockets socket)
{
	//get size of the incomming message
	message_t sizeMessage;
	fSockets[socket]->recv(&sizeMessage);
	int numberOfRecords;
	istringstream iss(static_cast<char*>(sizeMessage.data()));
	iss >> numberOfRecords;

	if(numberOfRecords==0){cout<<"MANAGER -- list is empty"<<endl;}
	
	fSockets[socket]->send(*(new message_t()));//receive empty message just to keep req-rep order
	
//get list of events
	message_t *response = new message_t(sizeof(serverListStruct)*numberOfRecords);
	fSockets[socket]->recv(response);
	
	vector<serverListStruct> receivedList(static_cast<serverListStruct*>(response->data()), static_cast<serverListStruct*>(response->data()) + numberOfRecords);

	if (response) {delete response;}
	return receivedList;
}

AliESDEvent* AliStorageEventManager::GetEvent(storageSockets socket,int timeout)
{
  message_t* message = new message_t();

  try
  {
    fSockets[socket]->recv(message);
  }
  catch (const zmq::error_t& e)
  {
    cout<<"MANAGER -- "<<e.what()<<endl;
    return NULL;
  }
	
  TBufferFile *mess = new TBufferFile(TBuffer::kRead,
                                      message->size()+sizeof(UInt_t),
                                      message->data());
  mess->InitMap();
  mess->ReadClass();// get first the class stored in message
  mess->SetBufferOffset(sizeof(UInt_t) + sizeof(kMESS_OBJECT));
  mess->ResetMap();
	
  AliESDEvent* data = (AliESDEvent*)(mess->ReadObjectAny(AliESDEvent::Class()));

  if (data)
  {
    data->GetStdContent();
    if(message){delete message;}
    return data; 
    /*
    TTree* tree= new TTree("esdTree", "esdTree");
    data->WriteToTree(tree);
    tree->Fill();
    AliESDEvent* event= new AliESDEvent();
    event->ReadFromTree(tree);
    tree->GetEntry(0);
    if(data){delete data;}
    if(tree){delete tree;}
    if(message){delete message;}
    return event;
    */
  }
  else
  {
    if(message){delete message;}
    return NULL;
  }
}

struct serverRequestStruct* AliStorageEventManager::GetServerStruct(storageSockets socket)
{
	struct serverRequestStruct *request = new struct serverRequestStruct;
	message_t *requestMessage = new message_t();
	fSockets[socket]->recv(requestMessage);
	request = static_cast<struct serverRequestStruct*>(requestMessage->data());
	return request;
}

struct clientRequestStruct* AliStorageEventManager::GetClientStruct(storageSockets socket)
{
	struct clientRequestStruct *request = new struct clientRequestStruct;
	message_t *requestMessage = new message_t();
	fSockets[socket]->recv(requestMessage);
	request = static_cast<struct clientRequestStruct*>(requestMessage->data());
	return request;
}

bool AliStorageEventManager::GetBool(storageSockets socket)
{
	message_t *response = new message_t();
	fSockets[socket]->recv(response);
	char *result = (char*)response->data();
	
	if(!strcmp("true",result)){return true;}
	else{return false;}
}

long AliStorageEventManager::GetLong(storageSockets socket)
{
	message_t *responseMessage = new message_t();
	fSockets[socket]->recv(responseMessage);
    
    long result = 0;
    
    if(responseMessage)
    {
        result = (long)atoi(static_cast<char*>(responseMessage->data()));
        delete responseMessage;
    }
	return result;
}

