#include "AliStorageEventManager.h"
#include "AliNetMessage.h"

#include <iostream>
#include <sstream>

#include <TList.h>
#include <TStreamerInfo.h>

#include "zmq.hpp"

using namespace zmq;
using namespace std;

AliStorageEventManager::AliStorageEventManager(){}
AliStorageEventManager::~AliStorageEventManager(){}

void __freeBuffer (void *data, void *hint)
{
    free(data);
}

void AliStorageEventManager::Send(AliESDEvent *event, socket_t *socket)
{
	AliNetMessage tmess(kMESS_OBJECT);
	tmess.Reset();
	tmess.WriteObject(event);

	int bufSize = tmess.BufferSize();
	char* buf = new char[bufSize];
	memcpy(buf, (char*)tmess.Buffer(), bufSize);

	message_t message(buf, bufSize, __freeBuffer, NULL);
	//fwrite(mess.Buffer(), sizeof(char), bufSize, stdout);
	
	socket->send(message);

	
	//publisher.Send(tmess);
			/*
	TMessage tmess(kMESS_OBJECT);
	tmess.Reset();
	tmess.WriteObject(event);
	TMessage::EnableSchemaEvolutionForAll(kTRUE);
	SendStreamerInfos(&tmess, socket);
	int bufsize = tmess.Length();
	char* buf = (char*) malloc(bufsize * sizeof(char));
	memcpy(buf, tmess.Buffer(), bufsize);
	zmq::message_t message((void*)buf, bufsize, 0, 0);
	socket->send(message);*/
}

void AliStorageEventManager::Send(vector<serverListStruct> list, socket_t *socket)
{
	//send size of the struct first
	int numberOfRecords = list.size();
	message_t message(20);
	snprintf ((char *)message.data(), 20 ,"%d",numberOfRecords);

	socket->send(message);
	if(numberOfRecords==0)return;
	socket->recv((new message_t));//empty message just to keep req-rep order

	//prepare message with event's list
	char *buffer = reinterpret_cast<char*> (&list[0]);
	message_t *reply = new message_t((void*)buffer,
			      sizeof(serverListStruct)*numberOfRecords,0);
	socket->send(*reply);

	if(reply){delete reply;}
}

void AliStorageEventManager::Send(struct serverRequestStruct *request,zmq::socket_t *socket)
{
	char *buffer = (char*)(request);
	message_t *requestMessage = new message_t((void*)buffer,
					   sizeof(struct serverRequestStruct)
					   +sizeof(struct listRequestStruct)
					   +sizeof(struct eventStruct),0);	
	socket->send(*requestMessage);
}

void AliStorageEventManager::Send(struct clientRequestStruct *request,zmq::socket_t *socket)
{
	char *buffer = (char*)(request);
	message_t *requestMessage = new message_t((void*)buffer,
					   sizeof(struct clientRequestStruct),0);	
	socket->send(*requestMessage);
}

void AliStorageEventManager::Send(long message,zmq::socket_t *socket)
{
	stringstream streamBuffer;
	streamBuffer << message;
	string stringBuffer = streamBuffer.str();
	char *buffer = (char*)stringBuffer.c_str();
	message_t *replyMessage = new message_t((void*)buffer,sizeof(long),0);
	socket->send(*replyMessage);
	delete replyMessage;
	streamBuffer.str(string());
	streamBuffer.clear();
}

void AliStorageEventManager::Send(bool message,zmq::socket_t *socket)
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
	message_t *replyMessage = new message_t((void*)buffer,sizeof(char*),0);
	socket->send(*replyMessage);
	delete replyMessage;
}

vector<serverListStruct> AliStorageEventManager::GetServerListVector(socket_t *socket)
{
	//get size of the incomming message
	message_t sizeMessage;
	socket->recv(&sizeMessage);
	int numberOfRecords;
	istringstream iss(static_cast<char*>(sizeMessage.data()));
	iss >> numberOfRecords;

	if(numberOfRecords==0){cout<<"MANAGER -- list is empty"<<endl;}
	
	socket->send(*(new message_t()));//receive empty message just to keep req-rep order
	
//get list of events
	message_t *response = new message_t(sizeof(serverListStruct)*numberOfRecords);
	socket->recv(response);
	
	vector<serverListStruct> receivedList(static_cast<serverListStruct*>(response->data()), static_cast<serverListStruct*>(response->data()) + numberOfRecords);

	return receivedList;
}

AliESDEvent* AliStorageEventManager::GetEvent(socket_t *socket)
{
	message_t *message = new message_t();
	
	socket->recv(message);
	int bufSize = (int)message->size();
		
	char* buf = new char[bufSize];
	memcpy(buf, (char*)message->data(), bufSize);
	
	AliNetMessage *mess = new AliNetMessage(buf, bufSize);
	
	/*
	message_t* message = RecvStreamerInfos(socket);
	int64_t more;
	size_t more_size = sizeof more;
    
	socket->getsockopt(ZMQ_RCVMORE, &more, &more_size );
	socket->recv(message);
	TBufferFile *mess = new TBufferFile(TBuffer::kRead,
					    message->size()+sizeof(UInt_t),
					    message->data());
	mess->InitMap();
	mess->ReadClass();// get first the class stored in message
	mess->SetBufferOffset(sizeof(UInt_t) + sizeof(kMESS_OBJECT));
	mess->ResetMap();
	

	//message_t *mess;
	//socket->recv(message);
	*/
	AliESDEvent* data = (AliESDEvent*)(mess->ReadObjectAny(AliESDEvent::Class()));
	if (data)
	{
		TTree* tree= new TTree("esdTree", "esdTree");
		data->WriteToTree(tree);
		tree-> Fill();
		AliESDEvent* event= new AliESDEvent();
		event->ReadFromTree(tree);
		tree->GetEntry(0);
		delete data;
		delete tree;
		if(message){delete message;}
		return event;    	    	
	}
	else
	{
		if(message){delete message;}
		return NULL;
	}
}

message_t* AliStorageEventManager::RecvStreamerInfos(socket_t *socket)
{
	message_t *message = new message_t;
	socket->recv(message);
	TBufferFile *mess = new TBufferFile(TBuffer::kRead,
					    message->size()+2*sizeof(UInt_t),
					    message->data());
	if(!mess)
	{
		if(message){delete message;}
		return NULL;
	}
	mess->InitMap();
	mess->ReadClass();     // get first the class stored in message
	mess->SetBufferOffset(sizeof(UInt_t) + sizeof(TMessage));
	mess->ResetMap();
	TList *list = (TList*)mess->ReadObjectAny(TList::Class());
	if(list==0 || list->IsEmpty())
	{
		return message;
	}
	if(message){delete message;}
	TIter next(list);
	TStreamerInfo *info;
	TObjLink *lnk = list->FirstLink();
	// First call BuildCheck for regular class
	while (lnk)
	{
		info = (TStreamerInfo*)lnk->GetObject();
		TObject *element = info->GetElements()->UncheckedAt(0);
		Bool_t isstl = element && strcmp("This",element->GetName())==0;
		if (!isstl)
		{
			info->BuildCheck();
			if (gDebug > 0)
			{
				Info("RecvStreamerInfos",
				     "importing TStreamerInfo: %s, version = %d",
				     info->GetName(),
				     info->GetClassVersion());
			}
		}
		lnk = lnk->Next();
	}
	lnk = list->FirstLink();
	while (lnk)
	{
		info = (TStreamerInfo*)lnk->GetObject();
		TObject *element = info->GetElements()->UncheckedAt(0);
		Bool_t isstl = element && strcmp("This",element->GetName())==0;
		if (isstl)
		{
			info->BuildCheck();
			if (gDebug > 0)
			{
				Info("RecvStreamerInfos",
				     "importing TStreamerInfo: %s, version = %d",
				     info->GetName(),
				     info->GetClassVersion());
			}
		}
		lnk = lnk->Next();
	}
	if(list){delete list;}
	if(mess){delete mess;}
	if(lnk){delete lnk;}
	return message;
}

void AliStorageEventManager::SendStreamerInfos(TMessage* mess, socket_t *socket)
{
	TList* infos = mess->GetStreamerInfos();
   
	TIter next(infos);
	TStreamerInfo *info;
	TList *minilist = 0;
	while ((info = (TStreamerInfo*)next()))
	{
		//Int_t uid = info->GetNumber();
		if (!minilist) minilist = new TList();
         
		minilist->Add(info);
	}
      
	if (minilist)
	{
		TMessage messinfo(kMESS_STREAMERINFO);
		messinfo.WriteObject(minilist);
		delete minilist;
		if (messinfo.GetStreamerInfos())
			messinfo.GetStreamerInfos()->Clear();
          
		int bufsize = messinfo.Length();
        	char* buf = (char*) malloc(bufsize * sizeof(char));
		memcpy(buf, messinfo.Buffer(), bufsize);

        	// send!
		zmq::message_t message((void*)buf, bufsize, 0, 0);
                     
		if (socket->send(message, ZMQ_SNDMORE))
		{
			Warning("SendStreamerInfos", "problems sending TStreamerInfo's ...");
		}
	}

	return;
}


