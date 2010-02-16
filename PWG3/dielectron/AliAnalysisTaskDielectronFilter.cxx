/*************************************************************************
* Copyright(c) 1998-2009, ALICE Experiment at CERN, All rights reserved. *
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

///////////////////////////////////////////////////////////////////////////
//                                                                       //
//                        Basic Analysis Task                            //
//                                                                       //
///////////////////////////////////////////////////////////////////////////

#include <TChain.h>

#include <AliLog.h>
#include <AliAODHandler.h>
#include <AliAnalysisManager.h>
#include <AliVEvent.h>

#include "AliDielectron.h"
#include "AliDielectronHistos.h"
#include "AliAnalysisTaskDielectronFilter.h"

ClassImp(AliAnalysisTaskDielectronFilter)

//_________________________________________________________________________________
AliAnalysisTaskDielectronFilter::AliAnalysisTaskDielectronFilter() :
  AliAnalysisTaskSE(),
  fDielectron(0)
{
  //
  // Constructor
  //
}

//_________________________________________________________________________________
AliAnalysisTaskDielectronFilter::AliAnalysisTaskDielectronFilter(const char *name) :
  AliAnalysisTaskSE(name),
  fDielectron(0)
{
  //
  // Constructor
  //
  DefineInput(0,TChain::Class());
  DefineOutput(1, THashList::Class());
}

//_________________________________________________________________________________
void AliAnalysisTaskDielectronFilter::Init()
{
  // Initialization
  if (fDebug > 1) AliInfo("Init() \n");
  
// require AOD handler
  AliAODHandler *aodH = (AliAODHandler*)((AliAnalysisManager::GetAnalysisManager())->GetOutputEventHandler());
  if (!aodH) Fatal("Init", "No AOD handler. Halting.");
  
//require dielectron framework 
  if (!fDielectron) {
    Error("Init","Dielectron framework class required. Please create and instance with proper cuts and set it via 'SetDielectron' before executing this task!!!");
    return;
  }
  
  aodH->AddFilteredAOD("AliAOD.Dielectron.root", "DielectronEvents");
//   AddAODBranch("AliDielectronCandidates",fDielectron->GetPairArraysPointer(),"deltaAOD.Dielectron.root");
}

//_________________________________________________________________________________
void AliAnalysisTaskDielectronFilter::UserExec(Option_t *)
{
  //
  // Main loop. Called for every event
  //

  if (!fDielectron) return;
  
  //bz for AliKF
  Double_t bz = InputEvent()->GetMagneticField();
  AliKFParticle::SetField( bz );
  
  fDielectron->Process(InputEvent());
  fDielectron->FillHistograms();

  if(fDielectron->HasCandidates()){
    AliAODExtension *extDielectron = dynamic_cast<AliAODHandler*>
      ((AliAnalysisManager::GetAnalysisManager())->GetOutputEventHandler())->GetFilteredAOD("AliAOD.Dielectron.root");
    extDielectron->SelectEvent();
    //see if dielectron candidate branch exists, if not create is
    TTree *t=extDielectron->GetTree();
    if (!t->GetBranch("dielectrons")){
      t->Bronch("dielectrons","TObjArray",fDielectron->GetPairArraysPointer());
    }
  }
  
  PostData(1, const_cast<THashList*>(fDielectron->GetHistogramList()));
}

