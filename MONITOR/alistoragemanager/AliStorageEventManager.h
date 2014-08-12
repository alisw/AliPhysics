#ifndef AliStorageEventManager_H
#define AliStorageEventManager_H

#include "AliESDEvent.h"
#include "AliStorageTypes.h"

#include <vector>
#include <string>
#include <sstream>

#include <TMessage.h>

namespace zmq
{
	class context_t;
	class socket_t;
}

class AliStorageEventManager
{
public:
	static AliStorageEventManager* GetEventManagerInstance();

	void Send(std::vector<serverListStruct> list,storageSockets socket);
	void Send(struct serverRequestStruct *request,storageSockets socket);
	bool Send(struct clientRequestStruct *request,storageSockets socket,int timeout = -1);
	void Send(AliESDEvent *event,storageSockets socket);
	void Send(long message,storageSockets socket);
	void Send(bool message,storageSockets socket);
	void SendAsXml(AliESDEvent *event,storageSockets socket);
	
	std::vector<serverListStruct> GetServerListVector(storageSockets socket);
	AliESDEvent* GetEvent(storageSockets socket,int timeout=-1,TTree **tmpTree=0);
	struct serverRequestStruct* GetServerStruct(storageSockets socket);
	struct clientRequestStruct* GetClientStruct(storageSockets socket);
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
	
	zmq::context_t *fContexts[7];
	zmq::socket_t *fSockets[7];

	AliStorageEventManager(const AliStorageEventManager&);
	AliStorageEventManager& operator=(const AliStorageEventManager&);
};

#endif
