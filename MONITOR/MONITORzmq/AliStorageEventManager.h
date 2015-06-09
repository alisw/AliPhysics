#ifndef AliStorageEventManager_H
#define AliStorageEventManager_H

#include "AliESDEvent.h"
#include "AliStorageTypes.h"

#include <vector>
#include <string>
#include <sstream>

#include <TMessage.h>
#include <TFile.h>

namespace zmq
{
	class context_t;
	class socket_t;
  class message_t;
}

class AliStorageEventManager
{
public:
	static AliStorageEventManager* GetEventManagerInstance();

	void Send(std::vector<serverListStruct> list,storageSockets socket);
	bool Send(struct serverRequestStruct *request,storageSockets socket,int timeout = -1);
	bool Send(struct clientRequestStruct *request,storageSockets socket,int timeout = -1);
	void Send(AliESDEvent *event,storageSockets socket);
	void Send(TFile *file,storageSockets socket);
	void Send(struct recPointsStruct *files,storageSockets socket);
	void Send(long message,storageSockets socket);
	void Send(bool message,storageSockets socket);
	void SendAsXml(AliESDEvent *event,storageSockets socket);
	
	std::vector<serverListStruct> GetServerListVector(storageSockets socket,int timeout=-1);
	AliESDEvent* GetEvent(storageSockets socket,int timeout=-1);
  zmq::message_t* GetMessage(storageSockets socket,int timeout=-1);
	TFile* GetFile(storageSockets socket,int timeout=-1);
	struct recPointsStruct* GetFiles(storageSockets socket,int timeout=-1);
	struct serverRequestStruct* GetServerStruct(storageSockets socket);
	struct clientRequestStruct* GetClientStruct(storageSockets socket,int timeout=-1);
	long GetLong(storageSockets socket);
	bool GetBool(storageSockets socket);


	bool CreateSocket(storageSockets socket);
private:
	AliStorageEventManager();
	~AliStorageEventManager();

	static AliStorageEventManager *fManagerInstance;

	std::string fStorageServer;
	std::string fEventServer;
	int fStorageServerPort;
	int fStorageClientPort;
	int fEventServerPort;
	int fXmlServerPort;
	int fItsPointsServerPort;
	
	zmq::context_t *fContexts[NUMBER_OF_SOCKETS];
	zmq::socket_t *fSockets[NUMBER_OF_SOCKETS];

	AliStorageEventManager(const AliStorageEventManager&);
	AliStorageEventManager& operator=(const AliStorageEventManager&);
};

#endif
