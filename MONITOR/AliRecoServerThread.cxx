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
  : AliThreadedSocket(context, AliThreadedSocket::kWRITE),
  	fReco(0),
    fCond(0)
{
	fReco = reco;
}

AliRecoServerThread::~AliRecoServerThread()
{
	Wait();
	Stop();
}

Bool_t AliRecoServerThread::Start(const char* endpoint)
{
	fHost = endpoint;
	
	return AliThreadedSocket::Start();
}

void AliRecoServerThread::RunThrdWrite()
{
	TThread::SetCancelAsynchronous();
	TThread::SetCancelOn();
	
 // generate a publish socket
	AliSocket publisher(fContext, ZMQ_PUB);
	publisher.Bind(fHost);
	
  if(fReco==0) return;
  
  AliESDEvent* event;
  
	fReco->Begin(NULL);
  if (fReco->GetAbort() != TSelector::kContinue) return;
  
  fReco->SlaveBegin(NULL);
	if (fReco->GetAbort() != TSelector::kContinue) return;
  
  //******* The loop over events
    Int_t iEvent = 0;
    while ( fReco->HasNextEventAfter(iEvent) ) {
      // check if process has enough resources 
      if (!fReco->HasEnoughResources(iEvent)) break;
      Bool_t status = fReco->ProcessEvent(iEvent);
      
      if (status)
      {
		    event = fReco->GetESDEvent();
		    
      	AliNetMessage tmess(kMESS_OBJECT);
  			tmess.Reset();
  			tmess.WriteObject(event);
  				
				publisher.Send(tmess);

   			sleep(1);
     }
      else {
        fReco->Abort("ProcessEvent",TSelector::kAbortFile);
      }
      		
      fReco->CleanProcessedEvent();
      /*
      if(fCond->TimedWaitRelative(500)==0){
				break;
				}			*/
      iEvent++;
    }
    fReco->SlaveTerminate();
    if (fReco->GetAbort() != TSelector::kContinue) return;
    fReco->Terminate();
    if (fReco->GetAbort() != TSelector::kContinue) return;
  
}
