#ifndef ALIANALYSISTASKFILTER_H
#define ALIANALYSISTASKFILTER_H
 
/* Copyright(c) 1998-1999, ALICE Experiment at CERN, All rights reserved. *
 * See cxx source for full Copyright notice                               */

/* $Id$ */

//////////////////////////////////////////////////////////////////////////
//
//  Base class for filtering friends
//
//////////////////////////////////////////////////////////////////////////

#include "AliAnalysisTask.h"
class AliVEvent;
class AliESDEvent;
class AliESDfriend;
class AliInputEventHandler;
class AliESDfriendTrack;

class TTree;

class AliAnalysisTaskFilter : public AliAnalysisTask
{
 public:
	AliAnalysisTaskFilter();
	AliAnalysisTaskFilter(const char* name);
	AliAnalysisTaskFilter(const AliAnalysisTaskFilter& obj);
	AliAnalysisTaskFilter& operator=(const AliAnalysisTaskFilter& other);
	virtual ~AliAnalysisTaskFilter() {;}

	// Implementation of interface methods
	virtual void ConnectInputData(Option_t *option = "");
	virtual void CreateOutputObjects();
	virtual void Exec(Option_t* option);
	virtual void SetDebugLevel(Int_t level) {fDebug = level;}
	virtual void Init() {;}
	virtual void UserCreateOutputObjects()  {;}
	virtual void UserExec(Option_t* /*option*/) {;}

	// To be implemented by user
	virtual Bool_t UserSelectESDfriendForCurrentEvent(){return kTRUE;}

	// Getters
	virtual Int_t         DebugLevel()  {return fDebug;     }
	virtual AliVEvent*    InputEvent()  {return fInputEvent;}
	virtual AliESDEvent*  ESDEvent()    {return fOutputESD; }
	virtual AliESDfriend* ESDfriend()   {return fOutputESDfriend; }
	virtual TTree*        OutputTree()  {return fTreeE;     }
	virtual Long64_t      Entry()       {return fEntry;     }
	virtual const char*   CurrentFileName();

	// To add a friend track
	void AddFriendTrackAt(AliESDfriendTrack* t, Int_t index);

 protected:
	Int_t                 fDebug;           //  Debug flag
	Int_t                 fEntry;           //  Current entry in the chain
	AliVEvent*            fInputEvent;      //! VEvent Input
	AliInputEventHandler* fInputHandler;    //! Input Handler
	AliESDEvent*          fOutputESD;       //! ESD out 
	AliESDfriend*         fOutputESDfriend; //! ESD friend out 
	TTree*                fTreeE;           //  ESD output Tree
	
	ClassDef(AliAnalysisTaskFilter, 1); // Analysis task for filtering friends
};

#endif
