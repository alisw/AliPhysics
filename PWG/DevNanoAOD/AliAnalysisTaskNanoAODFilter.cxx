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
  fTrkrep(0),
  fVarList(""),
  fVarListHead(""),
  fVarListHeader_fTC(""),
  fEvtCuts(0),
  fTrkCuts(0),
  fSetter(0),
  fSaveCutsFlag(0),
  fSaveAODZDC(kFALSE),
  fSaveVzero(kFALSE),
  fInputArrayName(""),
  fOutputArrayName("")
{
  // Dummy constructor ALWAYS needed for I/O.
}

//________________________________________________________________________
AliAnalysisTaskNanoAODFilter::AliAnalysisTaskNanoAODFilter(const char *name, Bool_t saveCutsFlag) // All data members should be initialised here
  :AliAnalysisTaskSE(name),
   fMCMode(0),
   fTrkrep(0),
   fVarList(""),
   fVarListHead(""),
   fVarListHeader_fTC(""),
   fEvtCuts(0),
   fTrkCuts(0),
   fSetter(0),
   fSaveCutsFlag(saveCutsFlag),
   fSaveAODZDC(kFALSE),
   fSaveVzero(kFALSE),
   fInputArrayName(""),
   fOutputArrayName("")

{
  // Constructor
  if(fSaveCutsFlag) {
    DefineOutput(1, AliAnalysisCuts::Class());
    DefineOutput(2, AliAnalysisCuts::Class());
  }


}

//________________________________________________________________________
AliAnalysisTaskNanoAODFilter::~AliAnalysisTaskNanoAODFilter()
{
  // Destructor. Clean-up the output list, but not the histograms that are put inside
  // (the list is owner and will clean-up these histograms). Protect in PROOF case.
}

//________________________________________________________________________
void AliAnalysisTaskNanoAODFilter::UserCreateOutputObjects()
{
  // Create histograms
  // Called once (on the worker node)
      
  if(fSaveCutsFlag) {
    PostData(1, fEvtCuts); 
    PostData(2, fTrkCuts); 
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
  

  AliNanoAODReplicator * rep = new AliNanoAODReplicator("NanoAODReplicator",
							"remove non interesting tracks, "
							"writes special tracks array tracks",
							fVarList,
							fVarListHead,
							fTrkCuts,
							fMCMode);

     
  cout<<"rep: "<<rep<<endl;
  rep->SetCustomSetter(fSetter);
  if (fSaveVzero) rep->SetVzero(1);
  if (fSaveAODZDC) rep->SetAODZDC(1);
  if (fVarListHeader_fTC) rep->SetVarListHeaderStringVariable(fVarListHeader_fTC);
  if (!fInputArrayName.IsNull()) rep->SetInputArrayName(fInputArrayName);
  if (!fOutputArrayName.IsNull()) rep->SetOutputArrayName(fOutputArrayName);

  std::cout << "SETTER: " << fSetter << " " << rep->GetCustomSetter() << std::endl;

  ext->DropUnspecifiedBranches(); // all branches not part of a FilterBranch call (below) will be dropped
      
  ext->FilterBranch("tracks",rep);
  ext->FilterBranch("vertices",rep);  
  ext->FilterBranch("header",rep);  
            
  if ( fMCMode > 0 ) 
    {
      // MC branches will be copied (if present), as they are, but only
      // for events with at least one muon. 
      // For events w/o muon, mcparticles array will be empty and mcheader will be dummy
      // (e.g. strlen(GetGeneratorName())==0)
      
      ext->FilterBranch("mcparticles",rep);
      ext->FilterBranch("mcHeader",rep);
    }
}

void AliAnalysisTaskNanoAODFilter::Init()
{

  // Initialization
  AddFilteredAOD("AliAOD.NanoAOD.root", "NanoAODTracksEvents");

}



//________________________________________________________________________
void AliAnalysisTaskNanoAODFilter::UserExec(Option_t *) 
{
  // Main loop
  Long64_t ientry = Entry();
  if(fDebug)printf("Nano AOD Filter: Analysing event # %5d\n", (Int_t) ientry);

  AliAODEvent *lAODevent=(AliAODEvent*)InputEvent();
  
    
  AliAODHandler* handler = dynamic_cast<AliAODHandler*>(AliAnalysisManager::GetAnalysisManager()->GetOutputEventHandler());

  if(fEvtCuts && !fEvtCuts->IsSelected(lAODevent)) return;// FIXME: should event cuts be called here or in the branch replicator? Do we get duplicated events if we skip here (arrays not reset in the branch replicator?)

  if ( handler ){
    AliAODExtension *extNanoAOD = handler->GetFilteredAOD("AliAOD.NanoAOD.root");
   if ( extNanoAOD ) {				
     extNanoAOD->SetEvent(lAODevent);
     extNanoAOD->SelectEvent();
     extNanoAOD->FinishEvent();
   }
  }


  if(fSaveCutsFlag) {
    PostData(1, fEvtCuts); 
    PostData(2, fTrkCuts); 
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
