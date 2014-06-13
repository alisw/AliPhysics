// Author: Mihai Niculescu 2013

/**************************************************************************
 * Copyright(c) 1998-2013, ALICE Experiment at CERN, all rights reserved. *
 * See http://aliceinfo.cern.ch/Offline/AliRoot/License.html for          *
 * full copyright notice.                                                 *
 **************************************************************************/
 
#ifndef __AliRecoServerThread_H__
#define __AliRecoServerThread_H__

#include <TQObject.h>
#include <RQ_OBJECT.h>
#include <TMutex.h>
#include <TCondition.h>

#include "AliThreadedSocket.h"

class TCondition;
class TThread;

class AliReconstruction;
class AliESDEvent;

class AliRecoServerThread : public AliThreadedSocket
{
public:
  AliRecoServerThread(zmq::context_t *context, AliReconstruction* reco);
  virtual ~AliRecoServerThread();


	Bool_t Start(const char* endpoint);

	const char* GetHost() const { return fHost.Data(); }	
	AliReconstruction*	GetReconstruction() { return fReco; }
	TCondition*					Condition() { return fCond; }
	
private:
	static void* RunThrdWrite(void* arg);
	
	AliReconstruction* fReco;

	// local	
  TCondition*						fCond; // condition whether to stop reco/clean exit thread
  TString fHost;

private:
  AliRecoServerThread(const AliRecoServerThread&);            // Not implemented
  AliRecoServerThread& operator=(const AliRecoServerThread&); // Not implemented
  
public:
  
  ClassDef(AliRecoServerThread, 0);  
};
#endif /* __AliReconstructionThread_H__ */

