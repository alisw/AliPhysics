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

using std::cout;
using std::endl;

ClassImp(AliAnalysisTaskNanoAODFilter)


//________________________________________________________________________
AliAnalysisTaskNanoAODFilter::AliAnalysisTaskNanoAODFilter() // All data members should be initialised here
:AliAnalysisTaskSE(),
  fMCMode(0),
  fReplicator(0),
  fEvtCuts(0),
  fTrkCuts(0),
  fV0Cuts(0),
  fSaveCutsFlag(0),
  fInputArrayName(""),
  fOutputArrayName("")
{
  // Dummy constructor ALWAYS needed for I/O.
}

//________________________________________________________________________
AliAnalysisTaskNanoAODFilter::AliAnalysisTaskNanoAODFilter(const char *name, Bool_t saveCutsFlag) // All data members should be initialised here
  :AliAnalysisTaskSE(name),
   fMCMode(0),
   fReplicator(0),
   fEvtCuts(0),
   fTrkCuts(0),
   fV0Cuts(0),
   fSaveCutsFlag(saveCutsFlag),
   fInputArrayName(""),
   fOutputArrayName("")

{
  fReplicator = new AliNanoAODReplicator("NanoAODReplicator", "remove non interesting tracks, writes special tracks array tracks");
  
  // Constructor
  if(fSaveCutsFlag) {
    DefineOutput(1, AliAnalysisCuts::Class());
    DefineOutput(2, AliAnalysisCuts::Class());
    DefineOutput(3, AliAnalysisCuts::Class());
  }
}

//________________________________________________________________________
AliAnalysisTaskNanoAODFilter::~AliAnalysisTaskNanoAODFilter()
{
  // Destructor. Clean-up the output list, but not the histograms that are put inside
  // (the list is owner and will clean-up these histograms). Protect in PROOF case.
  
  delete fReplicator;
}

//________________________________________________________________________
void AliAnalysisTaskNanoAODFilter::UserCreateOutputObjects()
{
  // Create histograms
  // Called once (on the worker node)
      
  if(fSaveCutsFlag) {
    if (fEvtCuts)
      PostData(1, fEvtCuts); 
    if (fTrkCuts)
      PostData(2, fTrkCuts); 
    if (fV0Cuts)
      PostData(3, fV0Cuts);
  }
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
  
  fReplicator->SetTrackCuts(fTrkCuts);
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
  
  if(fEvtCuts && !fEvtCuts->IsSelected(lAODevent)) return;// FIXME: should event cuts be called here or in the branch replicator? Do we get duplicated events if we skip here (arrays not reset in the branch replicator?)

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


//________________________________________________________________________
void AliAnalysisTaskNanoAODFilter::Terminate(Option_t *) 
{
  // print some debug info

}

void AliAnalysisTaskNanoAODFilter::FinishTaskOutput() {

  // We save here the user info //

  AliAODHandler* handler = dynamic_cast<AliAODHandler*>(AliAnalysisManager::GetAnalysisManager()->GetOutputEventHandler());
  AliAODExtension *extNanoAOD = handler->GetFilteredAOD("AliAOD.NanoAOD.root");
  extNanoAOD->GetTree()->GetUserInfo()->Add(AliNanoAODTrackMapping::GetInstance());

}
