/**************************************************************************
 * Copyright(c) 1998-2007, ALICE Experiment at CERN, All rights reserved. *
 *                                                                        *
 * Author: The ALICE Off-line Project.                                    *
 * Contributors are mentioned in the code where appropriate.              *
 *                                                                        *
 * Permission to use, copy, modify and distribute this software and its   *
 * documentation strictly for non-commercial purposes is hereby granted   *
 * without fee, provided that the above copyright notice appears in all   *
 * copies and that both the copyright notice and this permission notice   *
 * appear in the supporting documentation. The authors make no claims     *
 * about the suitability of this software for any purpose. It is          *
 * provided "as is" without express or implied warranty.                  *
 **************************************************************************/

//-------------------------------------------------------------------------
//     Event handler for reconstruction
//     Author: Andrei Gheata, CERN
//-------------------------------------------------------------------------

#include "AliHLTVEventInputHandler.h"
#include "AliVCuts.h"
#include "AliVEvent.h"
#include "TObjArray.h"
#include "TList.h"
#include "AliAnalysisTask.h"
#include "AliVfriendEvent.h"
#include "TTree.h"

ClassImp(AliHLTVEventInputHandler)

//______________________________________________________________________________
AliHLTVEventInputHandler::AliHLTVEventInputHandler() 
  : AliVEventHandler()
  , fEvent(NULL)
  , fFriendEvent(NULL)
  , fFakeTree(new TTree)
{
// default constructor
  fFakeTree->SetDirectory(0);
}

//______________________________________________________________________________
AliHLTVEventInputHandler::AliHLTVEventInputHandler(const char* name, const char* title) 
  : AliVEventHandler(name,title)
  , fEvent(NULL)
  , fFriendEvent(NULL)
  , fFakeTree(new TTree)
{
// Named constructor
  fFakeTree->SetDirectory(0);
}

//______________________________________________________________________________
AliHLTVEventInputHandler::AliHLTVEventInputHandler(AliHLTVEventInputHandler& that) 
  : AliVEventHandler(that)
  , fEvent(that.fEvent)
  , fFriendEvent(that.fFriendEvent)
  , fFakeTree(new TTree)
{
// dummy cpy constructor
  fFakeTree->SetDirectory(0);
}

AliHLTVEventInputHandler::~AliHLTVEventInputHandler()
{
    delete fFakeTree;
}

//______________________________________________________________________________
Bool_t AliHLTVEventInputHandler::Init(TTree* /*tree*/,  Option_t* /*opt*/)
{
// Initialisation necessary for each new tree. In reco case this is once.
  //Printf("----> AliHLTVEventInputHandler::Init"); 
  //Printf("<---- AliHLTVEventInputHandler::Init"); 

   return kTRUE;
}  
//______________________________________________________________________________
Bool_t AliHLTVEventInputHandler::BeginEvent(Long64_t)
{
// Called at the beginning of every event   

  //Printf("----> HLTTestInputHandler: BeginEvent: now fEvent is %p", fEvent);

  //Printf("----> HLTTestInputHandler: at the end of BeginEvent: now fEvent is %p", fEvent);
  return kTRUE;
}     

//______________________________________________________________________________
Bool_t AliHLTVEventInputHandler::AliHLTVEventInputHandler::FinishEvent()
{
  // Called at the end of every event   
  //when we are processing an AliESDEvent we have to detach the previously
  //attached friend as is is owned by us, not AliESDEvent
  //Printf("----> HLTTestInputHandler: FinishEvent");
  TList* list = fEvent->GetList();
  if (list) { list->Remove(fFriendEvent); }

  return kTRUE;
}     

//______________________________________________________________________________
Bool_t AliHLTVEventInputHandler::InitTaskInputData(AliVEvent* esdEvent, AliVfriendEvent* friendEvent, TObjArray* arrTasks) {

// Method to propagte to all the connected tasks the HLT event.
// The method gets the list of tasks from the manager

  //Printf("----> AliHLTVEventInputHandler::InitTaskInpuData: Setting the event...");
  SetEvent(esdEvent);
  SetVFriendEvent(friendEvent);
  esdEvent->SetFriendEvent(friendEvent);
  // set transient pointer to event inside tracks
  fEvent->ConnectTracks();
  //Printf("----> AliHLTVEventInputHandler::InitTaskInpuData: ...Event set: fEvent = %p; friend = %p", fEvent, friendEvent);
  for (Int_t i = 0; i < arrTasks->GetEntries(); i++){
    AliAnalysisTask* t = (AliAnalysisTask*)(arrTasks->At(i));
  //Printf("  ----> AliHLTVEventInputHandler: calling ConnectInputData() for task %i, %s",i,t->GetName());
    t->ConnectInputData("");
  }
  return kTRUE;
}

//______________________________________________________________________________
TTree* AliHLTVEventInputHandler::GetTree() const
{
  return fFakeTree;
}
