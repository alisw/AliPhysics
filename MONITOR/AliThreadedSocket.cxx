#include <zmq.hpp>

#include <TThread.h>

#include "AliNetMessage.h"
#include "AliSocket.h"
#include "AliThreadedSocket.h"

ClassImp(AliThreadedSocket)
AliThreadedSocket::AliThreadedSocket(zmq::context_t *context, EOpenMode mode)
	: TQObject(),
	fThread(0),
	fContext(context),
	fOpenMode(mode)
{

}

AliThreadedSocket::~AliThreadedSocket()
{
	Stop();
}

Bool_t AliThreadedSocket::Start()
{
	if(!fThread){
		if(fOpenMode==READ)
  		fThread = new TThread("AliThreadedSocket", (void(*) (void *) ) &RunThrdRead, (void*)  this );
  	else
  		fThread = new TThread("AliThreadedSocket", (void(*) (void *) ) &RunThrdWrite,(void*)  this );
  		
  	if(fThread->Run()==0){ 
	  	Emit("Started()");
  		return kTRUE; 
  	}
	}
	
	return kFALSE;
}

Bool_t AliThreadedSocket::Stop()
{
	Emit("Stopped()");
	return kTRUE;
}

Bool_t AliThreadedSocket::Kill()
{
	if(fThread){
		if(fThread->Kill()!=0) return kFALSE;
		fThread->Delete();
		fThread=0;
		
		Emit("Stopped()");
		return kTRUE;
	}
}

void AliThreadedSocket::Continue()
{
	
}

zmq::context_t* AliThreadedSocket::GetContext() const
{
	return fContext;
}

TThread* AliThreadedSocket::GetThread() const
{
	return fThread;
}

void AliThreadedSocket::Started()
{
	Emit("Started()");
}

void AliThreadedSocket::Stopped()
{
	Emit("Stopped()");
}

void* AliThreadedSocket::RunThrdRead(void* arg)
{
	AliNetMessage* mess=0;
	AliThreadedSocket* thsock = (AliThreadedSocket*)arg;
	zmq::context_t* context = thsock->GetContext();
	
	AliSocket sock(context, ZMQ_SUB);
		
	do{
		sock.Recv(mess);
	}
	while(mess==0);
	
	thsock->Stopped();
}

void* AliThreadedSocket::RunThrdWrite(void* arg)
{
	AliNetMessage* mess=0;
	AliThreadedSocket* thsock = (AliThreadedSocket*)arg;
	zmq::context_t* context = thsock->GetContext();
	
	AliSocket sock(context, ZMQ_PUB);
	
	do{
		sock.Send(*mess);
	}
	while(1);
	
	thsock->Emit("Stopped()");
}
