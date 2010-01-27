#ifndef ALIESDHANDLER_H
#define ALIESDHANDLER_H
/* Copyright(c) 1998-2007, ALICE Experiment at CERN, All rights reserved. *
 * See cxx source for full Copyright notice                               */

/* $Id$ */

//-------------------------------------------------------------------------
//
//     Implementation of the Event Handler Interface for ESD
//
//-------------------------------------------------------------------------

#include "AliVEventHandler.h"

class AliESDEvent;
class AliESDfriend;
class TFile;
class TTree;

class AliESDHandler : public AliVEventHandler {
    
 public:
	AliESDHandler();
	AliESDHandler(const char* name, const char* title);
	virtual ~AliESDHandler();
	virtual void SetOutputFileName(const char* fname){fFileName = fname;}
	virtual const char* GetOutputFileName(){return fFileName.Data();}
	virtual Bool_t Init(Option_t* option);
	virtual Bool_t Init(TTree* /*tree*/, Option_t* /*option*/)  {return kTRUE;}
	virtual Bool_t BeginEvent(Long64_t /*entry*/){return kTRUE;}
	virtual Bool_t Notify() {return AliVEventHandler::Notify(); };
	virtual Bool_t Notify(const char * /* path */) {return kTRUE;}
	virtual Bool_t FinishEvent();
	virtual Bool_t Terminate();
	virtual Bool_t TerminateIO();
	
	AliESDEvent* GetESD()  {return fESDEvent;}
	AliESDfriend* GetESDfriend()  {return fesdf;}
	virtual TTree* GetTree() const {return fTreeE;}
	void CreateTree(Int_t flag);
	void CreateFriends(Int_t flag);
	void FillTree();
	void AddESDtoTreeUserInfo();
	void SetInputTree(TTree* /*tree*/) {;}
 private:

	AliESDHandler(const AliESDHandler&);             // Not implemented
	AliESDHandler& operator=(const AliESDHandler&);  // Not implemented
	
	AliESDEvent* fESDEvent; //! Pointer to the ESD event
	AliESDfriend* fesdf;    //! Pointer to the ESD friend
	TTree* fTreeE;          //! tree for ESD persistency
	TFile* fFileE;          //! Output file
	TFile* fFileEF;         //! Output file for friends
	TString fFileName;      //  Output file name
	
    ClassDef(AliESDHandler, 1)
};
#endif
