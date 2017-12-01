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
class AliVCuts;

class AliESDHandler : public AliVEventHandler {
    
 public:
	AliESDHandler();
	AliESDHandler(const char* name, const char* title);
	virtual ~AliESDHandler();
	virtual void SetOutputFileName(const char* fname){fFileName = fname;}
	virtual const char* GetOutputFileName() const {return fFileName.Data();}
	virtual Bool_t Init(Option_t* option);
	virtual Bool_t Init(TTree* /*tree*/, Option_t* /*option*/)  {return kTRUE;}
	virtual Bool_t GetEntry() {return kTRUE;}
	virtual Bool_t BeginEvent(Long64_t /*entry*/){fIsEventSelectedForFriends = kFALSE; return kTRUE;}
	virtual Bool_t Notify() {return AliVEventHandler::Notify(); };
	virtual Bool_t Notify(const char * /* path */) {return kTRUE;}
	virtual Bool_t FinishEvent();
	virtual Bool_t Terminate();
	virtual Bool_t TerminateIO();
	
	AliESDfriend* GetESDfriend()  {return fesdf;}
	virtual TTree* GetTree() const {return fTreeEF;}
	void FillTree();
	void SetInputTree(TTree* /*tree*/) {;}
	void SelectEventForFriends() {fIsEventSelectedForFriends = kTRUE;}
  virtual AliVCuts*    GetEventSelection() const {return NULL;}

 private:

	AliESDHandler(const AliESDHandler&);             // Not implemented
	AliESDHandler& operator=(const AliESDHandler&);  // Not implemented
	
	AliESDfriend* fesdf;    //! Pointer to the ESD friend
	TTree* fTreeEF;         //! Output tree for friends
	TFile* fFileEF;         //! Output file for friends
	TString fFileName;      //! Output file name for friends
	Bool_t fIsEventSelectedForFriends; //! flag to indicate if the event was selected to have the friends kept 

    ClassDef(AliESDHandler, 3)
};
#endif
