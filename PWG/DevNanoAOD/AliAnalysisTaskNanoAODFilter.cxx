/**************************************************************************
 * Copyright(c) 1998-1999, ALICE Experiment at CERN, All rights reserved. *
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

/* $Id$ */

// Generic task to produce nanoAOD
// Author: Michele Floris, michele.floris@cern.ch

#include "AliAnalysisTaskNanoAODFilter.h"

#include "Riostream.h"
#include "TChain.h"
#include "TTree.h"
#include "TH1F.h"
#include "TH2F.h"
#include "TCanvas.h"
#include "TList.h"

#include "AliAnalysisTaskSE.h"
#include "AliAnalysisManager.h"
#include "AliStack.h"
#include "AliESDtrackCuts.h"
#include "AliESDEvent.h"
#include "AliESDInputHandler.h"
#include "AliAODEvent.h"
#include "AliAODTrack.h"
#include "AliMCEvent.h"
#include "AliESDtrack.h"
#include "AliAODHandler.h"
#include "AliNanoAODReplicator.h"
#include "AliNanoAODTrackMapping.h"
#include "AliNanoAODTrack.h"
#include "AliNanoFilterNormalisation.h"
#include "AliMultSelectionTask.h"

ClassImp(AliAnalysisTaskNanoAODFilter)


//________________________________________________________________________
AliAnalysisTaskNanoAODFilter::AliAnalysisTaskNanoAODFilter() // All data members should be initialised here
:AliAnalysisTaskSE(),
  fUseAliEventCuts(false),
  fEventCuts(),
  fMCMode(0),
  fReplicator(0),
  fEvtCuts(),
  fQAOutput(0),
  fSaveCutsFlag(kFALSE),
  fInputArrayName(""),
  fOutputArrayName(""),
  fNormalisation(0x0)

{
  // Dummy constructor ALWAYS needed for I/O.
}

//________________________________________________________________________
AliAnalysisTaskNanoAODFilter::AliAnalysisTaskNanoAODFilter(const char *name, Bool_t saveCutsFlag) // All data members should be initialised here
  :AliAnalysisTaskSE(name),
   fUseAliEventCuts(false),
   fEventCuts(),
   fMCMode(0),
   fReplicator(0),
   fEvtCuts(0),
   fQAOutput(0),
   fSaveCutsFlag(saveCutsFlag),
   fInputArrayName(""),
   fOutputArrayName(""),
   fNormalisation(0x0)

{
  // Constructor

  fReplicator = new AliNanoAODReplicator("NanoAODReplicator", "remove non interesting tracks, writes special tracks array tracks");
  fQAOutput = new TList();
  fQAOutput->SetOwner(kTRUE); 
  
  DefineOutput(1, AliNanoFilterNormalisation::Class());
  DefineOutput(2, TList::Class());
}

//________________________________________________________________________
AliAnalysisTaskNanoAODFilter::~AliAnalysisTaskNanoAODFilter()
{
  // Destructor. Clean-up the output list, but not the histograms that are put inside
  // (the list is owner and will clean-up these histograms). Protect in PROOF case.
  
  delete fReplicator;
  if (fQAOutput) 
    delete fQAOutput;
  delete fNormalisation;
}

//________________________________________________________________________
void AliAnalysisTaskNanoAODFilter::UserCreateOutputObjects()
{
  // Create histograms
  // Called once (on the worker node)
      
  if (fSaveCutsFlag) {
    for (std::list<AliAnalysisCuts*>::iterator it = fEvtCuts.begin(); it != fEvtCuts.end(); ++it)
      fQAOutput->Add(*it);
    PostData(2, fQAOutput);
  }
  
  std::string normName = std::string(fName) + "_scaler";
  fNormalisation = new AliNanoFilterNormalisation(normName.data(), normName.data());
  PostData(1, fNormalisation);
}

void AliAnalysisTaskNanoAODFilter::AddFilteredAOD(const char* aodfilename, const char* title)
{
  // The replicator is added to the extension

  AliAODHandler *aodH = (AliAODHandler*)((AliAnalysisManager::GetAnalysisManager())->GetOutputEventHandler());
  if (!aodH) AliFatal("No AOD handler");
  AliAODExtension* ext = aodH->AddFilteredAOD(aodfilename,title);
    
  if (!ext){
    AliFatal("Cannot get extension");
  }
  
  fReplicator->SetMCMode(fMCMode);
     
  if (!fInputArrayName.IsNull()) fReplicator->SetInputArrayName(fInputArrayName);
  if (!fOutputArrayName.IsNull()) fReplicator->SetOutputArrayName(fOutputArrayName);

  ext->DropUnspecifiedBranches(); // all branches not part of a FilterBranch call (below) will be dropped
      
  ext->FilterBranch("tracks",fReplicator);
  ext->FilterBranch("vertices",fReplicator);  
  ext->FilterBranch("header",fReplicator);  
            
  if ( fMCMode > 0 ) 
    {
      // MC branches will be copied (if present), as they are, but only
      // for events with at least one muon. 
      // For events w/o muon, mcparticles array will be empty and mcheader will be dummy
      // (e.g. strlen(GetGeneratorName())==0)
      
      ext->FilterBranch("mcparticles",fReplicator);
      ext->FilterBranch("mcHeader",fReplicator);
    }
}

void AliAnalysisTaskNanoAODFilter::Init()
{
  // Initialization
  AddFilteredAOD("AliAOD.NanoAOD.root", "NanoAODTracksEvents");
}

void AliAnalysisTaskNanoAODFilter::UserExec(Option_t *) 
{
  // Main loop
  Long64_t ientry = Entry();
  if(fDebug)printf("Nano AOD Filter: Analysing event # %5d\n", (Int_t) ientry);

  AliAODEvent *lAODevent= dynamic_cast<AliAODEvent*> (InputEvent());
  if (!lAODevent)
    lAODevent = AODEvent(); // On the fly ESD filtering
  if (!lAODevent)
    AliFatal("No input event");

  fNormalisation->FillCandidate(kTRUE, kTRUE, kTRUE, kTRUE, 0);
  
  for (std::list<AliAnalysisCuts*>::iterator it = fEvtCuts.begin(); it != fEvtCuts.end(); ++it)
    if (!((*it)->IsSelected(lAODevent))) 
      return;

  if (fUseAliEventCuts) {
    fEventCuts.AcceptEvent(lAODevent);
    double mult = AliMultSelectionTask::IsINELgtZERO(lAODevent) ? fEventCuts.GetCentrality() : -0.5;
    fNormalisation->FillSelected(fEventCuts.CheckNormalisationMask(AliEventCuts::kTriggeredEvent),
                                 fEventCuts.CheckNormalisationMask(AliEventCuts::kPassesNonVertexRelatedSelections),
                                 fEventCuts.CheckNormalisationMask(AliEventCuts::kHasReconstructedVertex),
                                 fEventCuts.CheckNormalisationMask(AliEventCuts::kPassesAllCuts),
                                 mult);
  } else
    fNormalisation->FillSelected(kTRUE, kTRUE, kTRUE, kTRUE, 0);

  AliAODHandler* handler = dynamic_cast<AliAODHandler*>(AliAnalysisManager::GetAnalysisManager()->GetOutputEventHandler());
  if ( handler ){
    AliAODExtension *extNanoAOD = handler->GetFilteredAOD("AliAOD.NanoAOD.root");
   if ( extNanoAOD ) {				
     extNanoAOD->SetEvent(lAODevent);
     extNanoAOD->SelectEvent();
     extNanoAOD->FinishEvent();
   }
  }
}

void AliAnalysisTaskNanoAODFilter::Terminate(Option_t *) 
{
}

void AliAnalysisTaskNanoAODFilter::FinishTaskOutput() {

  // We save here the user info

  AliAODHandler* handler = dynamic_cast<AliAODHandler*>(AliAnalysisManager::GetAnalysisManager()->GetOutputEventHandler());
  AliAODExtension *extNanoAOD = handler->GetFilteredAOD("AliAOD.NanoAOD.root");

  // copy production version info
  AliVEventHandler* inputHandler = AliAnalysisManager::GetAnalysisManager()->GetInputEventHandler();
  if (inputHandler->GetUserInfo()) {
    TObject* prodInfo = inputHandler->GetUserInfo()->FindObject("alirootVersion");
    if (prodInfo)
      extNanoAOD->GetTree()->GetUserInfo()->Add(prodInfo->Clone());
  }
  
  Printf("****************************************************************");
  extNanoAOD->GetTree()->GetUserInfo()->Add(AliNanoAODTrackMapping::GetInstance(fReplicator->GetVarListTrack()));
  AliNanoAODTrackMapping::GetInstance()->Print();
  Printf("****************************************************************");
  
  extNanoAOD->GetTree()->GetUserInfo()->Add(fNormalisation->Clone());
}

void AliAnalysisTaskNanoAODFilter::AddPIDField(AliNanoAODTrack::ENanoPIDResponse response, AliPID::EParticleType particle)
{
  TString list(fReplicator->GetVarListTrack());
  list += ",";
  list += AliNanoAODTrack::GetPIDVarName(response, particle);
  fReplicator->SetVarListTrack(list);
}
