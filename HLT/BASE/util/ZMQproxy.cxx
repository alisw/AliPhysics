#include "zmq.h"
#include <iostream>
#include "AliHLTDataTypes.h"
#include "AliHLTComponent.h"
#include "AliHLTMessage.h"
#include "TClass.h"
#include "TMap.h"
#include "TPRegexp.h"
#include "TObjString.h"
#include "TList.h"
#include "TMessage.h"
#include "TRint.h"
#include "TApplication.h"
#include <time.h>
#include <unistd.h>
#include <string>
#include <map>
#include "AliZMQhelpers.h"

//this is meant to become a class, hence the structure with global vars etc.
//Also the code is rather flat - it is a bit of a playground to test ideas.
//TODO structure this at some point, e.g. introduce a SIMPLE unified way of handling
//zmq payloads, maybe a AliZMQmessage class which would by default be multipart and provide
//easy access to payloads based on topic or so (a la HLT GetFirstInputObject() etc...)

//methods
int ProcessOptionString(TString arguments);
int InitZMQ();
void* work(void* param);
int Run();

//configuration vars
TString fZMQconfigIN   = "";
TString fZMQconfigOUT  = "";
TString fZMQconfigMON  = "";

Bool_t  fSendOnMerge = kTRUE;
Bool_t  fResetOnSend = kFALSE;

//ZMQ stuff
void* fZMQcontext = NULL;    //ze zmq context

void* fZMQmon = NULL;        //the request-reply socket, here we request the merged data
void* fZMQout = NULL;        //the monitoring socket, here we publish a copy of the data
void* fZMQin  = NULL;        //the in socket - entry point for the data to be merged.
int inType = -1;
int outType = -1;
int monType = -1;

int fZMQtimeout = -1;
int fZMQmaxQueueSize = 100;
int fSleep = 0;
bool fVerbose = false;

const char* fUSAGE =
  "ZMQproxy: a simple monitored ZMQ proxy\n"
  "caveat: using a REQ socket causes a custom request to be sent;\n"
  "        only the reply is forwardedto the backend.\n"
  "        For request forwarding use DEALER-ROUTER.\n"
  "options:\n"
  " -in : socket in\n"
  " -out : socket out\n"
  " -mon : monitor socket\n"
  " -sleep : sleep between polls (ms) (if in is a REQ socket)\n"
  " -timeout : timeout for a poll (ms)\n"
  ;

int capture(void* capture, zmq_msg_t* msg, int more = 0)
{
  int rc = 0;
  if (!capture) return 0;
  zmq_msg_t ctrl;
  rc = zmq_msg_init(&ctrl);
  if (rc<0) return -1;
  rc = zmq_msg_copy(&ctrl, msg);
  if (rc<0) {
    zmq_msg_close(&ctrl);
    return -1;
  }
  rc = zmq_msg_send(&ctrl, capture, more?ZMQ_SNDMORE:0);
  if (rc<0) {
    zmq_msg_close(&ctrl);
  }
  return rc;
}

int forward(void* from, void* to, void* mon, zmq_msg_t* msg)
{
  int rc = 0;
  int more = 0;
  size_t moresize = sizeof(more);
  while (true)
  {
    rc = zmq_msg_recv(msg,from,0);
    if (rc<0) return -1;
    rc = zmq_getsockopt(from, ZMQ_RCVMORE, &more, &moresize);
    if (from!=mon) {
      rc = capture(mon,msg,more);
    }
    if (rc<0) return -1;
    rc = zmq_msg_send(msg,to,more?ZMQ_SNDMORE:0);
    if (rc<0) return -1;
    if (more==0) break;
  }
  if (fVerbose) printf("forwarded...\n");
  return 0;
}

//_______________________________________________________________________________________
int Run()
{
  int rc = 0;

  bool interrunpted = false;
  while (!interrunpted)
  {
    //send the requests first
    //it is OK to block if there is no remote end
    //the use case is to go from REQ->PUB for the event display
    //otherwise the reconnect can take a while
    if (inType==ZMQ_REQ) {
      if (fVerbose) printf("requesting on %p %s\n",fZMQin, fZMQconfigIN.Data());
      rc = alizmq_msg_send("","",fZMQin,0);
    }
    if (outType==ZMQ_REQ) {
      if (fVerbose) printf("requesting on %s\n",fZMQconfigOUT.Data());
      alizmq_msg_send("","",fZMQout,0);
    }

    //poll socket readiness for sending
    Int_t noutSockets=2;
    zmq_pollitem_t outSockets[] = { 
      { fZMQin, 0, ZMQ_POLLOUT, 0 },
      { fZMQout, 0, ZMQ_POLLOUT, 0 },
      { fZMQmon, 0, ZMQ_POLLOUT, 0 },
    };
    rc = zmq_poll(outSockets, noutSockets, fZMQtimeout); //poll outSockets
    if (rc==-1 && errno==ETERM)
    {
      //this can only happen if the context was terminated, one of the inSockets are
      //not valid or operation was interrupted
      Printf("zmq_poll (out) was interrupted! rc = %i, %s", rc, zmq_strerror(errno));
      break;
    }
    //in
    if (outSockets[0].revents & ZMQ_POLLOUT)
    {
      if (fVerbose) printf("socket (%p) %s signals ZMQ_POLLOUT\n",fZMQin, fZMQconfigIN.Data());
    }
    //out
    if (outSockets[1].revents & ZMQ_POLLOUT)
    {
      if (fVerbose) printf("socket %p (%s) signals ZMQ_POLLOUT\n", fZMQout, fZMQconfigOUT.Data());
    }
    //mon
    if (outSockets[2].revents & ZMQ_POLLOUT)
    {
      if (fVerbose) printf("socket %p (%s) signals ZMQ_POLLOUT\n", fZMQmon, fZMQconfigMON.Data());
    }
    
    //poll incoming data
    Int_t ninSockets=3;
    zmq_pollitem_t inSockets[] = { 
      { fZMQin, 0, ZMQ_POLLIN, 0 },
      { fZMQout, 0, ZMQ_POLLIN, 0 },
      { fZMQmon, 0, ZMQ_POLLIN, 0 },
    };

    if (fVerbose) printf("starting poll in\n");
    int pollrc = zmq_poll(inSockets, ninSockets, fZMQtimeout); //poll inSockets
    if (pollrc==-1 && errno==ETERM)
    {
      //this can only happen if the context was terminated, one of the inSockets are
      //not valid or operation was interrupted
      Printf("zmq_poll (in) was interrupted! rc = %i, %s", rc, zmq_strerror(errno));
      break;
    }

    zmq_msg_t msg;
    zmq_msg_init(&msg);
    //in
    if (inSockets[0].revents & ZMQ_POLLIN)
    {
      if (fVerbose) printf("data in on %s\n",fZMQconfigIN.Data());
      if (outSockets[1].events & ZMQ_POLLOUT)
      {
        rc = forward(fZMQin,fZMQout,fZMQmon,&msg);
      }
    }
    //out
    if (inSockets[1].revents & ZMQ_POLLIN)
    {
      if (fVerbose) printf("data in on %s\n",fZMQconfigOUT.Data());
      if (outSockets[0].events & ZMQ_POLLOUT)
      {
        rc = forward(fZMQout,fZMQin,fZMQmon,&msg);
      }
    }
    //mon
    if (inSockets[2].revents & ZMQ_POLLIN)
    {
      if (fVerbose) printf("data in on %s\n",fZMQconfigMON.Data());
      if (outSockets[2].events & ZMQ_POLLOUT)
      {
        rc = forward(fZMQmon,fZMQmon,NULL,&msg);
      }
    }
    zmq_msg_close(&msg);
    
    //if we time out (waiting for a response) reinit the REQ socket(s)
    if (pollrc==0)
    {
      if (inType==ZMQ_REQ) {
        if (fVerbose) printf("no reply from %s in %i ms, server died?\n",
            fZMQconfigIN.Data(), fZMQtimeout);
        rc = alizmq_socket_init(fZMQin, fZMQcontext, fZMQconfigIN.Data(), 
            -1, fZMQmaxQueueSize);
        if (fVerbose) printf("rc of reinit %i\n",rc);
      }
      if (outType==ZMQ_REQ) {
        if (fVerbose) printf("no reply from %s in %i ms, server died?\n",
            fZMQconfigOUT.Data(), fZMQtimeout);
        rc = alizmq_socket_init(fZMQout, fZMQcontext, fZMQconfigOUT.Data(), 
            -1, fZMQmaxQueueSize);
        if (fVerbose) printf("rc of reinit %i\n",rc);
      }
      if (monType==ZMQ_REQ) {
        if (fVerbose) printf("no reply from %s in %i ms, server died?\n",
            fZMQconfigMON.Data(), fZMQtimeout);
        rc = alizmq_socket_init(fZMQmon, fZMQcontext, fZMQconfigMON.Data(), 
            -1, fZMQmaxQueueSize);
        if (fVerbose) printf("rc of reinit %i\n",rc);
      }
    }

    if (fSleep>0) usleep(fSleep);
  }

  return rc;
}

//_______________________________________________________________________________________
int InitZMQ()
{
  //init or reinit stuff
  int rc = 0;
  inType = alizmq_socket_init(fZMQin,  fZMQcontext, fZMQconfigIN.Data(), -1, fZMQmaxQueueSize);
  outType = alizmq_socket_init(fZMQout, fZMQcontext, fZMQconfigOUT.Data(), -1, fZMQmaxQueueSize);
  monType = alizmq_socket_init(fZMQmon, fZMQcontext, fZMQconfigMON.Data(), -1, fZMQmaxQueueSize);
  printf("socket init (in,out,mon): %i, %i, %i, %s, %s, %s\n",inType,outType,monType,alizmq_socket_name(inType),alizmq_socket_name(outType),alizmq_socket_name(monType));
  return rc;
}

//______________________________________________________________________________
int ProcessOptionString(TString arguments)
{
  //process passed options
  aliStringVec* options = AliOptionParser::TokenizeOptionString(arguments);
  int nOptions = 0;
  for (aliStringVec::iterator i=options->begin(); i!=options->end(); ++i)
  {
    const TString& option = i->first;
    const TString& value = i->second;
    if (option.EqualTo("ZMQconfigIN") || option.EqualTo("in"))
    {
      fZMQconfigIN = value;
    }
    else if (option.EqualTo("ZMQconfigOUT") || option.EqualTo("out"))
    {
      fZMQconfigOUT = value;
    }
    else if (option.EqualTo("ZMQconfigMON") || option.EqualTo("mon"))
    {
      fZMQconfigMON = value;
    }
    else if (option.EqualTo("Verbose"))
    {
      fVerbose = true;
    }
    else if (option.EqualTo("sleep"))
    {
      fSleep = value.Atoi() * 1000;
    }
    else if (option.EqualTo("timeout"))
    {
      fZMQtimeout = value.Atoi();
    }
    else
    {
      nOptions=-1;
      break;
    }
    nOptions++;
  }
  delete options; //tidy up

  return nOptions; 
}

//_______________________________________________________________________________________
int main(int argc, char** argv)
{
  int mainReturnCode=0;

  //process args
  TString argString = AliOptionParser::GetFullArgString(argc,argv);
  if (ProcessOptionString(argString)<=0)
  {
    printf("%s",fUSAGE);
    return 1;
  }

  //the context
  fZMQcontext = alizmq_context();

  //init stuff
  if (InitZMQ()<0) {
    Printf("failed init");
    return 1;
  }

  Run();

  //destroy ZMQ sockets
  zmq_close(fZMQmon);
  zmq_close(fZMQin);
  //zmq_close(fZMQout);
  zmq_ctx_destroy(fZMQcontext);
  return mainReturnCode;
}

