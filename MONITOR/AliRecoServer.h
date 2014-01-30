// Author: Mihai Niculescu 2013

/**************************************************************************
 * Copyright(c) 1998-2013, ALICE Experiment at CERN, all rights reserved. *
 * See http://aliceinfo.cern.ch/Offline/AliRoot/License.html for          *
 * full copyright notice.                                                 *
 **************************************************************************/
 
#ifndef __AliRecoServer_H__
#define __AliRecoServer_H__

#include <zmq.hpp>
#include <TObjString.h>
#include <TQObject.h>
#include <RQ_OBJECT.h>

class TEnv;
class AliCDBManager;
class AliReconstruction;
class AliRecoServerThread;

class AliRecoServer : public TQObject
{
	RQ_OBJECT("AliRecoServer")
public:
  AliRecoServer();
  virtual ~AliRecoServer();
  
  // Closes the server. The server will no longer listen/serve
  void Close();
  
  Bool_t IsListenning() const;
  
  Int_t GetRunId() const {  return fCurrentRunId; }
  const char* GetError() const;//TODO: not implemented
  
  Bool_t StartReconstruction(Int_t run, const char* input="mem://@*:");
  void  StopReconstruction();

  void ThreadFinished(Int_t status); // *SIGNAL*
  
private:
  Int_t RetrieveGRP(UInt_t run, TString &gdc);
  void SetupReco(const char* input);

  AliRecoServer(const AliRecoServer&);            // Not implemented
  AliRecoServer& operator=(const AliRecoServer&); // Not implemented

private:
	// thread shared
	zmq::context_t* fContext;
  AliReconstruction* fReco;
  AliCDBManager* fCDBman;

	// not shared
  Int_t										 fCurrentRunId;
  Bool_t fIsListenning;
  TEnv*								 fSettings;
  AliRecoServerThread* fRecoTh;

public:

  ClassDef(AliRecoServer, 0);
};

#endif /* __AliRecoServer_H__ */
