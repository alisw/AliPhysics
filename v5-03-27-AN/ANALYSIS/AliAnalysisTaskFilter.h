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
	virtual AliESDfriend* InputFriend() {return fInputESDfriend;}
	virtual AliESDfriend* ESDfriend()   {return fOutputESDfriend; }
	virtual TTree*        OutputTree()  {return fTreeEF;     }
	virtual Long64_t      Entry()       {return fEntry;     }
	virtual const char*   CurrentFileName();

	// To add/skip a friend track
	void AddFriendTrackAt(AliESDfriendTrack* t, Int_t index);
	void SkipFriendTrackAt(Int_t index);

 protected:
	Int_t                 fDebug;           //  Debug flag
	Int_t                 fEntry;           //  Current entry in the chain
	AliVEvent*            fInputEvent;      //! VEvent Input
	AliInputEventHandler* fInputHandler;    //! Input Handler
	AliESDfriend*         fOutputESDfriend; //! ESD friend out 
	TTree*                fTreeEF;          //  ESD friend output Tree
	AliESDfriend*         fInputESDfriend;  //! ESD friend input
	
	ClassDef(AliAnalysisTaskFilter, 3); // Analysis task for filtering friends
};

#endif
