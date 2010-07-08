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
#include <AliInputEventHandler.h>
#include <AliESDInputHandler.h>

#include "AliDielectron.h"
#include "AliDielectronMC.h"
#include "AliDielectronHistos.h"
#include "AliDielectronVarManager.h"
#include "AliAnalysisTaskDielectronFilter.h"

ClassImp(AliAnalysisTaskDielectronFilter)

//_________________________________________________________________________________
AliAnalysisTaskDielectronFilter::AliAnalysisTaskDielectronFilter() :
    AliAnalysisTaskSE(),
                      fDielectron(0),
                                  fSelectPhysics(kTRUE)
{
  //
  // Constructor
  //
}

//_________________________________________________________________________________
AliAnalysisTaskDielectronFilter::AliAnalysisTaskDielectronFilter(const char *name) :
    AliAnalysisTaskSE(name),
                      fDielectron(0),
                                  fSelectPhysics(kTRUE)
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
  fDielectron->Init();
  
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
  
  AliAnalysisManager *man=AliAnalysisManager::GetAnalysisManager();
  AliESDInputHandler *esdHandler=0x0;
  if ( (esdHandler=dynamic_cast<AliESDInputHandler*>(man->GetInputEventHandler())) && esdHandler->GetESDpid() ){
    AliDielectronVarManager::SetESDpid(esdHandler->GetESDpid());
  } else {
    //load esd pid bethe bloch parameters depending on the existance of the MC handler
    // yes: MC parameters
    // no:  data parameters
    if (!AliDielectronVarManager::GetESDpid()){
      if (AliDielectronMC::Instance()->HasMC()) {
        AliDielectronVarManager::InitESDpid();
      } else {
        AliDielectronVarManager::InitESDpid(1);
      }
    }
  }
  // Was event selected ?
  AliInputEventHandler* inputHandler = (AliInputEventHandler*) (man->GetInputEventHandler());
  Bool_t isSelected = kTRUE;
  if( fSelectPhysics && inputHandler && inputHandler->GetEventSelection() ) {
    isSelected = inputHandler->IsEventSelected();
  }
  
  if (!isSelected) return;
  
  //bz for AliKF
  Double_t bz = InputEvent()->GetMagneticField();
  AliKFParticle::SetField( bz );
  
  fDielectron->Process(InputEvent());
  
  if(fDielectron->HasCandidates()){
    //If input event is an AliESDevent
    // replace the references of the legs with the AOD references
    if(man->GetInputEventHandler()->IsA()==AliESDInputHandler::Class()){
      AliAODEvent *aod = ((AliAODHandler*)((AliAnalysisManager::GetAnalysisManager())->GetOutputEventHandler()))->GetAOD();
      TObjArray *obj = 0x0;
      AliAODTrack *leg1 = 0x0;
      AliAODTrack *leg2 = 0x0;
      for(Int_t i=0; i < 10; i++ ){
        obj = (TObjArray*)(*(fDielectron->GetPairArraysPointer()))->UncheckedAt(i);
        if(!obj) continue;
        for(int j=0;j<obj->GetEntriesFast();j++)
        {
          AliDielectronPair *pairObj = (AliDielectronPair*)obj->UncheckedAt(j);
          Int_t id1 = ((AliESDtrack*)pairObj->GetFirstDaughter())->GetID();
          Int_t id2 = ((AliESDtrack*)pairObj->GetSecondDaughter())->GetID();
          for(Int_t it=0;it<aod->GetNumberOfTracks();it++){
            if(aod->GetTrack(it)->GetID() == id1) leg1 = aod->GetTrack(it);
            if(aod->GetTrack(it)->GetID() == id2) leg2 = aod->GetTrack(it);
          }
          if(!leg1 || !leg2) continue;
          pairObj->SetRefFirstDaughter(leg1);
          pairObj->SetRefSecondDaughter(leg2);
        }
      }
    }
      
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

