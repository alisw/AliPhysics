// Author: Mihai Niculescu 2013

/**************************************************************************
 * Copyright(c) 1998-2013, ALICE Experiment at CERN, all rights reserved. *
 * See http://aliceinfo.cern.ch/Offline/AliRoot/License.html for          *
 * full copyright notice.                                                 *
 **************************************************************************/
 
#ifndef __AliRecoServer_H__
#define __AliRecoServer_H__

#include <TObjString.h>
#include <TQObject.h>
#include <RQ_OBJECT.h>
#include <TThread.h>

class TEnv;
class AliCDBManager;
class AliReconstruction;
class AliRecoServerThread;

class AliEventServerReconstruction : public TQObject
{
public:
	AliEventServerReconstruction();
	virtual ~AliEventServerReconstruction();

	Bool_t StartReconstruction(Int_t run, const char* input="mem://@*:");
	void  StopReconstruction();  
	
	// Closes the server. The server will no longer listen/serve
	void Close();
  
	Bool_t IsListenning() const{return fIsListenning;}
	Int_t GetRunId() const {return fCurrentRunId;}  
private:
    static void* Dispatch(void *arg){static_cast<AliEventServerReconstruction*>(arg)->ReconstructionHandle();}
	void ReconstructionHandle();

	Int_t RetrieveGRP(UInt_t run, TString &gdc);
	void SetupReco(const char* input);

	// thread shared
	AliReconstruction *fAliReco;
	AliCDBManager *fCDBmanager;
	Int_t fCurrentRunId;
	Bool_t fIsListenning;
	TEnv *fSettings;
	TString fHost;
	TThread *fRecoThread;

	AliEventServerReconstruction(const AliEventServerReconstruction&);
	AliEventServerReconstruction& operator=(const AliEventServerReconstruction&);
};
#endif
