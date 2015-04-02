#ifndef AliZMQManager_H
#define AliZMQManager_H

#include "AliESDEvent.h"
#include "AliStorageTypes.h"

#include <vector>
#include <string>
#include <sstream>

#include "zmq.h"

#include <TMessage.h>
#include <TFile.h>

class AliZMQManager
{
public:
	static AliZMQManager* GetInstance();
    void CreateSocket(storageSockets socket);
    void RecreateSocket(storageSockets socket);
    
	bool Send(std::vector<serverListStruct> list,storageSockets socket);
	bool Send(struct serverRequestStruct request,storageSockets socket);
	bool Send(struct clientRequestStruct *request,storageSockets socket);
	bool Send(AliESDEvent *event,storageSockets socket);
	bool Send(long message,storageSockets socket);
	bool Send(bool message,storageSockets socket);
    bool SendAsXml(AliESDEvent *event,storageSockets socket);
	
    bool Get(std::vector<serverListStruct>* &result,storageSockets socket);
    bool Get(AliESDEvent *result, storageSockets socket);
	bool Get(struct serverRequestStruct* &result, storageSockets socket);
	bool Get(struct clientRequestStruct *result, storageSockets socket);
	bool Get(long *result, storageSockets socket);
	bool Get(bool *result, storageSockets socket);

private:
	AliZMQManager();
	~AliZMQManager();
	static AliZMQManager *fManagerInstance;
    
    // ZMQ methods wrappers:
    bool zmqInit(zmq_msg_t *msg,size_t size=0);             // these methids return true on success
    bool zmqSend(zmq_msg_t *msg,void *socket,int flags);
    bool zmqRecv(zmq_msg_t *msg,void *socket,int flags);
    
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

	AliZMQManager(const AliZMQManager&);
	AliZMQManager& operator=(const AliZMQManager&);
};

#endif
