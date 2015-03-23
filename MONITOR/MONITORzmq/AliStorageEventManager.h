#ifndef AliStorageEventManager_H
#define AliStorageEventManager_H

#include "AliESDEvent.h"
#include "AliStorageTypes.h"

#include <vector>
#include <string>
#include <sstream>

#include "zmq.h"

#include <TMessage.h>
#include <TFile.h>

class AliStorageEventManager
{
public:
	static AliStorageEventManager* GetEventManagerInstance();
	bool CreateSocket(storageSockets socket);
    
	void Send(std::vector<serverListStruct> list,storageSockets socket);
	bool Send(struct serverRequestStruct *request,storageSockets socket,int timeout = -1);
	bool Send(struct clientRequestStruct *request,storageSockets socket,int timeout = -1);
	void Send(AliESDEvent *event,storageSockets socket);
	void Send(long message,storageSockets socket);
	void Send(bool message,storageSockets socket);
	void SendAsXml(AliESDEvent *event,storageSockets socket);
	
	std::vector<serverListStruct> GetServerListVector(storageSockets socket,int timeout=-1);
	AliESDEvent* GetEvent(storageSockets socket,int timeout=-1);
	struct serverRequestStruct* GetServerStruct(storageSockets socket);
	struct clientRequestStruct* GetClientStruct(storageSockets socket,int timeout=-1);
	long GetLong(storageSockets socket);
	bool GetBool(storageSockets socket);

private:
	AliStorageEventManager();
	~AliStorageEventManager();
	static AliStorageEventManager *fManagerInstance;
    
    // ZMQ methods wrappers:
    void zmqInit(zmq_msg_t *msg,size_t size=-1);
    void zmqSend(zmq_msg_t *msg,void *socket,int flags);
    void zmqRecv(zmq_msg_t *msg,void *socket,int flags);
    bool fZmqError; // if something went wrong in the above methods, this will be set to true
    
    // hostnames and ports read from config file:
	std::string fStorageServer;
	std::string fEventServer;
	int fStorageServerPort;
	int fStorageClientPort;
	int fEventServerPort;
	int fXmlServerPort;
	int fItsPointsServerPort;
	
    // ZMQ sockets and contexts:
	void *fContexts[NUMBER_OF_SOCKETS];
	void *fSockets[NUMBER_OF_SOCKETS];

	AliStorageEventManager(const AliStorageEventManager&);
	AliStorageEventManager& operator=(const AliStorageEventManager&);
};

#endif
