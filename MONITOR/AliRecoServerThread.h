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

class TCondition;
class TMessage;
class TThread;

class AliReconstruction;
class AliESDEvent;

namespace zmq{
	class context_t;
	class socket_t;
}

class AliRecoServerThread : public TQObject
{
public:
  AliRecoServerThread(zmq::context_t *context, AliReconstruction* reco);
  virtual ~AliRecoServerThread();

	Bool_t Start(const char* host);
	Int_t Stop();
	Bool_t ForceStop(); // imediate kill it, use it with rarely and with caution

	zmq::context_t*			GetContext() { return fContext; }
	AliReconstruction*	GetReconstruction() { return fReco; }
	const char*					GetHost() { return fHost.Data(); }
	TCondition*							Condition() { return fCond; }
	
	void Finished(Int_t status); // *SIGNAL*

private:
	static void* RunThreaded(void* arg);
	static void SendStreamerInfos(TMessage* mess, zmq::socket_t *sock);
	static void SendEvent(AliESDEvent* event, zmq::socket_t* socket);
	
	// shared
	zmq::context_t* fContext;
	AliReconstruction* fReco;

	// local	
	TString										fHost;
  TThread*								fThread;
  TCondition*						fCond; // condition whether to stop reco/clean exit thread

private:
  AliRecoServerThread(const AliRecoServerThread&);            // Not implemented
  AliRecoServerThread& operator=(const AliRecoServerThread&); // Not implemented
  
public:
  
  ClassDef(AliRecoServerThread, 0);  
};
#endif /* __AliReconstructionThread_H__ */

