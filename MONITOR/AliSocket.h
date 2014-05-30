#ifndef ALISOCKET_H
#define ALISOCKET_H

#include <TObject.h>
#include <TString.h>

class AliNetMessage;

namespace zmq
{
	class context_t;
	class socket_t;
}

class AliSocket : public TObject
{
public:
	AliSocket(zmq::context_t* context, int type);
	virtual ~AliSocket();
	
	void Bind(const char* endpoint);
	void Connect(const char* endpoint);
	void Subscribe(const char* filter);
	
  bool Recv(AliNetMessage *&mess, int flags = 0);
  bool Send(AliNetMessage &message, int flags = 0);
  
private:
	AliSocket(const AliSocket &); // Not implemented
	void operator=(const AliSocket &); // Not implemented
	
	zmq::context_t *fContext; //! the zmq context
	zmq::socket_t *fSocket; //! the socket
	TString fEndPoint; //!
	
  ClassDef(AliSocket, 0);  
};

#endif
