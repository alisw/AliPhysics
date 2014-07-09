#ifndef AliStorageEventManager_H
#define AliStorageEventManager_H

#include "AliESDEvent.h"
#include "AliStorageTypes.h"

#include <vector>

#include <TMessage.h>

namespace zmq
{
	class socket_t;
	class message_t;
}

class AliStorageEventManager
{
public:
	AliStorageEventManager();
	~AliStorageEventManager();

	void Send(std::vector<serverListStruct> list,zmq::socket_t *socket);
	void Send(struct serverRequestStruct *request,zmq::socket_t *socket);
	void Send(struct clientRequestStruct *request,zmq::socket_t *socket);
	void Send(AliESDEvent *event, zmq::socket_t *socket);
	void Send(long message,zmq::socket_t *socket);
	void Send(bool message,zmq::socket_t *socket);
	
	std::vector<serverListStruct> GetServerListVector(zmq::socket_t *socket);
	AliESDEvent* GetEvent(zmq::socket_t *socket);
	
private:
	void SendStreamerInfos(TMessage *mess, zmq::socket_t *socket);
	zmq::message_t* RecvStreamerInfos(zmq::socket_t *socket);
};

#endif
