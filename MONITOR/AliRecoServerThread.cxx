// Author: Mihai Niculescu 2013

/**************************************************************************
 * Copyright(c) 1998-2013, ALICE Experiment at CERN, all rights reserved. *
 * See http://aliceinfo.cern.ch/Offline/AliRoot/License.html for          *
 * full copyright notice.                                                 *
 **************************************************************************/

#include <RVersion.h>
#include <stdlib.h>

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

#include <AliNetMessage.h>
#include <AliSocket.h>

#include <zmq.hpp>

#include "AliRecoServerThread.h"

ClassImp(AliRecoServerThread)
AliRecoServerThread::AliRecoServerThread(zmq::context_t *context, AliReconstruction* reco)
  : AliThreadedSocket(context, AliThreadedSocket::WRITE),
  	fReco(0),
    fCond(0)
{
	fReco = reco;
}

AliRecoServerThread::~AliRecoServerThread()
{
	Stop();
}

Bool_t AliRecoServerThread::Start(const char* endpoint)
{
	fHost = endpoint;
	
	return AliThreadedSocket::Start();
}

void* AliRecoServerThread::RunThrdWrite(void* arg)
{
	TThread::SetCancelAsynchronous();
	TThread::SetCancelOn();
	
	AliRecoServerThread* recoTh = (AliRecoServerThread*)arg;
	
	const char* host = recoTh->GetHost();
	zmq::context_t* context = recoTh->GetContext();
	AliReconstruction* reco = recoTh->GetReconstruction();

 // generate a publish socket
	AliSocket publisher(context, ZMQ_PUB);
	publisher.Bind(host);
	
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
		    
      	AliNetMessage tmess(kMESS_OBJECT);
  			tmess.Reset();
  			tmess.WriteObject(event);
  				
				publisher.Send(tmess);

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
