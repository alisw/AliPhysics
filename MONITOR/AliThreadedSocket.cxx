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
	Wait();
	Stop();
}

Bool_t AliThreadedSocket::Start()
{
	if(!fThread){
  		fThread = new TThread("AliThreadedSocket", (void(*) (void *) ) &Dispatch, (void*)  this );
  
  	if(fThread->Run()==0){ 
	  	Emit("Started()");
  		return kTRUE; 
  	}
	}
	
	return kFALSE;
}

Bool_t AliThreadedSocket::Stop()
{
	if(fThread){
		fThread->Delete();
		fThread=0;
	}

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

void AliThreadedSocket::Wait()
{
	if(fThread && fThread->GetState()==TThread::kRunningState)
	{
		fThread->Join();
	}
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

void AliThreadedSocket::RunThrdRead()
{
	AliNetMessage* mess=0;
	AliSocket sock(fContext, ZMQ_SUB);
		
	do{
		sock.Recv(mess);
	}
	while(mess==0);
	
	Stopped();
}

void AliThreadedSocket::RunThrdWrite()
{
	AliNetMessage* mess=0;
	AliSocket sock(fContext, ZMQ_PUB);
	
	do{
		sock.Send(*mess);
	}
	while(1);
	
	Stopped();
}
