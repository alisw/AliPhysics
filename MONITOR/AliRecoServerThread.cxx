// Author: Mihai Niculescu 2013

/**************************************************************************
 * Copyright(c) 1998-2013, ALICE Experiment at CERN, all rights reserved. *
 * See http://aliceinfo.cern.ch/Offline/AliRoot/License.html for          *
 * full copyright notice.                                                 *
 **************************************************************************/

#include <RVersion.h>
#include <stdlib.h>

#include <zmq.hpp>

#include <TCondition.h>
#include <TBufferFile.h>
#include <TMessage.h>
#include <TObjArray.h>
#include <TStreamerInfo.h>
#include <TThread.h>


#include <AliESDEvent.h>
#include <AliESDfriend.h>
#include <AliRawReader.h>
#include <AliRunLoader.h>
#include <AliReconstruction.h>

#include "AliRecoServerThread.h"

ClassImp(AliRecoServerThread);
AliRecoServerThread::AliRecoServerThread(zmq::context_t *context, AliReconstruction* reco)
  : TQObject(),
		fContext(0),
  	fReco(0),
  	fHost("tcp://*:5051"),
    fThread(0),
    fCond(0)
{
	fContext = context;
	fReco = reco;
}

AliRecoServerThread::~AliRecoServerThread()
{
	Stop();
}

Bool_t AliRecoServerThread::Start(const char* host)
{
	if(!fThread){
  	fHost = host;
  	fCond = new TCondition(0);
  	fThread = new TThread("AliRecoServerThread", (void(*) (void *) ) &RunThreaded, (void*)  this );
  	fThread->Run();
  
  	return kTRUE;
	}
	
	return kFALSE;	
}

Int_t AliRecoServerThread::Stop()
{
	fCond->Signal();
 
  return 0;
}

Bool_t AliRecoServerThread::ForceStop()
{
	if(fThread){
		fThread->Kill();
		fThread->Delete();
		fThread=0;
		
		return kTRUE;
	}
	
	return kFALSE;
}

void AliRecoServerThread::Finished(Int_t status)
{
  Emit("Finished(Int_t)", status);
}

void AliRecoServerThread::SendStreamerInfos(TMessage* mess, zmq::socket_t *sock)
{
	//printf("Sending Streamer Infos....\n");

	// Check if TStreamerInfo must be sent. The list of TStreamerInfo of classes
   // in the object in the message is in the fInfos list of the message.
   // We send only the TStreamerInfos not yet sent on this socket.
	TList* infos = mess->GetStreamerInfos();
   
      TIter next(infos);
      TStreamerInfo *info;
      TList *minilist = 0;
      while ((info = (TStreamerInfo*)next())) {
         Int_t uid = info->GetNumber();
         if (!minilist) minilist = new TList();
         
         minilist->Add(info);
      }
      
      if (minilist) {
         TMessage messinfo(kMESS_STREAMERINFO);
         messinfo.WriteObject(minilist);
         delete minilist;
         if (messinfo.GetStreamerInfos())
            messinfo.GetStreamerInfos()->Clear();
          
					int bufsize = messinfo.Length();
        	char* buf = (char*) malloc(bufsize * sizeof(char));
          memcpy(buf, messinfo.Buffer(), bufsize);

        	// send!
          zmq::message_t message((void*)buf, bufsize, 0, 0);
                     
         if (sock->send(message, ZMQ_SNDMORE))
            Warning("SendStreamerInfos", "problems sending TStreamerInfo's ...");
      }

   return;
}

void AliRecoServerThread::SendEvent(AliESDEvent* event, zmq::socket_t* socket)
{
  if(!event) return;

  TMessage tmess(kMESS_OBJECT);
  tmess.Reset();
  tmess.WriteObject(event);

  TMessage::EnableSchemaEvolutionForAll(kTRUE);
  SendStreamerInfos(&tmess, socket);

  int bufsize = tmess.Length();
  char* buf = (char*) malloc(bufsize * sizeof(char));
  memcpy(buf, tmess.Buffer(), bufsize);

  // send!
  zmq::message_t message((void*)buf, bufsize, 0, 0);
  socket->send(message);

}


void* AliRecoServerThread::RunThreaded(void* arg)
{
	TThread::SetCancelAsynchronous();
	TThread::SetCancelOn();
	
	AliRecoServerThread* recoTh = (AliRecoServerThread*)arg;
	
	const char* host = recoTh->GetHost();
	zmq::context_t* context = recoTh->GetContext();
	AliReconstruction* reco   = recoTh->GetReconstruction();

	zmq::socket_t publisher(*context, ZMQ_PUB);
	publisher.bind(host);
	
  if(reco==0) return 0;
  
  AliESDEvent* event;
  
	reco->Begin(NULL);
  if (reco->GetAbort() != TSelector::kContinue) return 0;
  
  reco->SlaveBegin(NULL);
	if (reco->GetAbort() != TSelector::kContinue) return 0;
  
  //******* The loop over events
    Int_t iEvent = 0;
    while ( reco->HasNextEventAfter(iEvent) ) {
      // check if process has enough resources 
      if (!reco->HasEnoughResources(iEvent)) break;
      Bool_t status = reco->ProcessEvent(iEvent);
      
      if (status)
      {
				event = reco->GetESDEvent();
				SendEvent(event, &publisher);

   			sleep(1);
     }
      else {
        reco->Abort("ProcessEvent",TSelector::kAbortFile);
      }
      		
      reco->CleanProcessedEvent();
      if(recoTh->Condition()->TimedWaitRelative(500)==0){
				break;
			}			
      iEvent++;
    }
    reco->SlaveTerminate();
    if (reco->GetAbort() != TSelector::kContinue) return 0;
    reco->Terminate();
    if (reco->GetAbort() != TSelector::kContinue) return 0;
  
}
