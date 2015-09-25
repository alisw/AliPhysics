/// \class AliZMQManager
/// \brief ZMQ communication manager
///
/// This class handles communication based on ZMQ to provide easy to use interface.
/// Manager is capable of sending and receiving different data types and can be extended
/// for other types if needed.
///
/// Usage:
/// 1) Get manager's instance:
///     AliZMQManager *manager = AliZMQManager::GetInstance();
/// 2) Create a socket:
///     manager->CreateSocket(NAME_OF_SOCKET);
///
///     Names of avaliable sockets are stored in AliStorageTypes.h
/// 3) Use Send() and Get():
///
///    manager->Send(true,NAME_OF_SOCKET);
///    manager->Get(&longResponse,NAME_OF_SOCKET);
///
/// Communication's logic:
///     Every send and get call has a timeout of 5 seconds. Every Send(...) and Get(...) will
///     return bool, which indicates status (1-success,0-failed). It is your responsibility
///     to invent a communication logic that will handle all cases. Especially keep in mind,
///     that REQ-REP communication is sensitive to the state of socket, so if you successfully
///     sent massage on REQ socket and then you failed to receive a response, you cannot send
///     new request, as socket is expecting you to try to receive message again.
///
///     If you are sure you want to send new message after receive failed, you can always recreate
///     socket (it will be closed and created again in the default state). To do that, call:
///
///     manager->RecreateSocket(NAME_OF_SOCKET);
///
/// Adding new sockets:
///     If you need to add new socket, you need to do two things:
///     1) Add name of socket in AliStorageTypes.h (enum storageSockets),
///     2) Extend CreateSocket(...) method by adding another case in the switch (to be consistent, copy
///         one of existing cases and modify it).
///
/// \author Jeremi Niedziela <jeremi.niedziela@cern.ch>, Warsaw University of Technology, CERN

#ifndef AliZMQManager_H
#define AliZMQManager_H

#include "AliESDEvent.h"
#include "AliStorageTypes.h"

#include <vector>
#include <string>
#include <sstream>
#include <set>

#include <TMessage.h>
#include <TFile.h>

struct zmq_msg_t;

class AliZMQManager
{
public:
	static AliZMQManager* GetInstance();
    void Close();
    void CreateSocket(storageSockets socket);
    void RecreateSocket(storageSockets socket);
    
    // methods to send data (1 = success, 0 = timeout, -1 = socket closed):
	int Send(std::vector<serverListStruct> list,storageSockets socket);
    int Send(std::vector<string100> list,storageSockets socket);
	int Send(struct serverRequestStruct *request,storageSockets socket);
	int Send(struct clientRequestStruct *request,storageSockets socket);
	int Send(AliESDEvent *event,storageSockets socket);
	int Send(long message,storageSockets socket);
	int Send(bool message,storageSockets socket);
    int SendAsXml(AliESDEvent *event,storageSockets socket);
	
    // methods to receive data (1 = success, 0 = timeout, -1 = socket closed):
    int Get(std::vector<serverListStruct>* &result,storageSockets socket);
    int Get(std::vector<string100>* &result,storageSockets socket);
    int Get(AliESDEvent* &result, storageSockets socket);
	int Get(struct serverRequestStruct* &result, storageSockets socket);
	int Get(struct clientRequestStruct* &result, storageSockets socket);
	int Get(long *result, storageSockets socket);
	int Get(bool *result, storageSockets socket);

private:
	AliZMQManager();
	~AliZMQManager();
	static AliZMQManager *fManagerInstance; ///< single instance of AliZMQManager
    
    // ZMQ methods wrappers (1 = success, 0 = timeout, -1 = socket closed):
    bool zmqInit(zmq_msg_t *msg,size_t size=0);
    int zmqSend(zmq_msg_t *msg,void *socket,int flags);
    int zmqRecv(zmq_msg_t *msg,void *socket,int flags);
    
    // hostnames and ports read from config file:
	std::string fStorageServer; ///< hostname for Storage Manager
	std::string fEventServer;   ///< hostname for Online Reconstruction
	int fStorageServerPort;     ///< SM server thread's port
	int fStorageClientPort;     ///< SM client thread's port
	int fEventServerPort;       ///< Online Reconstruction port
	int fXmlServerPort;         ///< port for xml publisher
	
    // ZMQ sockets and contexts:
	void *fContext; ///< single ZMQ context for all sockets
	void *fSockets[NUMBER_OF_SOCKETS];  ///< array of all ZMQ sockets

    // copy constructor and assignment operator:
	AliZMQManager(const AliZMQManager&);
	AliZMQManager& operator=(const AliZMQManager&);
};

#endif
