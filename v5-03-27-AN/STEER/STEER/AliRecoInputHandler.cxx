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

#include "AliRecoInputHandler.h"
#include "AliVCuts.h"

ClassImp(AliRecoInputHandler)

//______________________________________________________________________________
AliRecoInputHandler::AliRecoInputHandler(const char* name, const char* title) 
  : AliESDInputHandler(name,title)
{
// Named constructor
}

//______________________________________________________________________________
Bool_t AliRecoInputHandler::Init(TTree* tree,  Option_t* opt)
{
// Initialisation necessary for each new tree. In reco case this is once.
   fAnalysisType = opt;
   fTree = tree;
   if (!fTree) return kFALSE;
   fNEvents = fTree->GetEntries();
   return kTRUE;
}  
//______________________________________________________________________________
Bool_t AliRecoInputHandler::BeginEvent(Long64_t)
{
// Called at the beginning of every event   
  static Bool_t called = kFALSE;
  if (!called && fEventCuts && IsUserCallSelectionMask())
     AliInfo(Form("The ESD input handler expects that the first task calls AliESDInputHandler::CheckSelectionMask() %s", fEventCuts->ClassName()));
  fNewEvent = kTRUE;
  called = kTRUE;
  return kTRUE;
}     
