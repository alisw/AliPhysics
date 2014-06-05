#include <cstring>

#include <zmq.hpp>

#include "AliNetMessage.h"
#include "AliSocket.h"

void __freeBuffer (void *data, void *hint)
{
    free(data);
}


ClassImp(AliSocket);
AliSocket::AliSocket(zmq::context_t *context,int type)
	: TObject(),
	fContext(context)
{
	fSocket = new zmq::socket_t(*context,type);
}

AliSocket::~AliSocket()
{
	fSocket->close();
}

void AliSocket::Bind(const char* endpoint)
{
	fEndPoint = endpoint;
	
	fSocket->bind(endpoint);
}

void AliSocket::Connect(const char* endpoint)
{
	fEndPoint = endpoint;
	
	fSocket->connect(endpoint);
}

void AliSocket::Subscribe(const char* filter)
{
	fSocket->setsockopt(ZMQ_SUBSCRIBE, filter, strlen(filter) );
}

bool AliSocket::Recv(AliNetMessage *&mess, int flags)
{
	zmq::message_t message;
	  
	if(!fSocket->recv(&message, flags))
		return false;	
	
	int bufSize = (int)message.size();
		
	// buffer will be adopted by AliNetMessage, no need to free it
	char* buf = new char[bufSize];
	memcpy(buf, (char*)message.data(), bufSize);
	
	mess = new AliNetMessage(buf, bufSize);

	return true;
}

bool AliSocket::Send(AliNetMessage &mess, int flags)
{
	//NOTE: this is already done by AliNetMessage, should we do this too?
	// send length of the message
	int bufSize = mess.BufferSize();
	
	// we need to copy it elsewhere because zmq takes ownership of the buffer data
	char* buf = new char[bufSize];
	memcpy(buf, (char*)mess.Buffer(), bufSize);

	zmq::message_t message(buf, bufSize, __freeBuffer, NULL);
		
	//fwrite(mess.Buffer(), sizeof(char), bufSize, stdout);
	
	return fSocket->send(message, flags);
}
