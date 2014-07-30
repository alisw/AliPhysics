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

#include "AliHLTTestInputHandler.h"
#include "AliVCuts.h"
#include "AliVEvent.h"
#include "TObjArray.h"
#include "AliAnalysisTask.h"

ClassImp(AliHLTTestInputHandler)

//______________________________________________________________________________
AliHLTTestInputHandler::AliHLTTestInputHandler(const char* name, const char* title) 
  : AliVEventHandler(name,title)
{
// Named constructor
}

//______________________________________________________________________________
Bool_t AliHLTTestInputHandler::Init(TTree* tree,  Option_t* opt)
{
// Initialisation necessary for each new tree. In reco case this is once.
  Printf("----> AliHLTTestInputHandler::Init"); 
  Printf("<---- AliHLTTestInputHandler::Init"); 

   return kTRUE;
}  
//______________________________________________________________________________
Bool_t AliHLTTestInputHandler::BeginEvent(Long64_t)
{
// Called at the beginning of every event   

  Printf("----> HLTTestInputHandler: BeginEvent: now fEvent is %p", fEvent);

  Printf("----> HLTTestInputHandler: at the end of BeginEvent: now fEvent is %p", fEvent);
  return kTRUE;
}     

//______________________________________________________________________________
Bool_t AliHLTTestInputHandler::InitTaskInputData(AliVEvent* esdEvent, AliESDfriend* friendEvent, TObjArray* arrTasks) {

// Method to propagte to all the connected tasks the HLT event.
// The method gets the list of tasks from the manager

  Printf("----> AliHLTTestInputHandler::InitTaskInpuData: Setting the event...");
  SetEvent(esdEvent);
  SetFriendEvent(friendEvent);
  // set transient pointer to event inside tracks
  fEvent->ConnectTracks();
  Printf("----> AliHLTTestInputHandler::InitTaskInpuData: ...Event set");
  for (Int_t i = 0; i < arrTasks->GetEntries(); i++){
    AliAnalysisTask* t = (AliAnalysisTask*)(arrTasks->At(i));
    t->ConnectInputData("");
  }
  return kTRUE;
}
