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

/* $Id$ */

///////////////////////////////////////////////////////////////////////////
//                                                                       //
//                        Basic Analysis Task                            //
//                                                                       //
///////////////////////////////////////////////////////////////////////////

#include <TChain.h>
#include <TH1D.h>

#include <AliLog.h>
#include <AliAODHandler.h>
#include <AliAODInputHandler.h>
#include <AliAnalysisManager.h>
#include <AliVEvent.h>
#include <AliInputEventHandler.h>
#include <AliESDInputHandler.h>
#include <AliAODInputHandler.h>

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
fSelectPhysics(kTRUE),
fTriggerMask(AliVEvent::kMB),
fEventStat(0x0),
fStoreLikeSign(kFALSE)
{
  //
  // Constructor
  //
}

//_________________________________________________________________________________
AliAnalysisTaskDielectronFilter::AliAnalysisTaskDielectronFilter(const char *name) :
AliAnalysisTaskSE(name),
fDielectron(0),
fSelectPhysics(kTRUE),
fTriggerMask(AliVEvent::kMB),
fEventStat(0x0),
fStoreLikeSign(kFALSE)
{
  //
  // Constructor
  //
  DefineInput(0,TChain::Class());
  DefineOutput(1, THashList::Class());
  DefineOutput(2, TH1D::Class());
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
void AliAnalysisTaskDielectronFilter::UserCreateOutputObjects()
{
  //
  // Initilise histograms
  //
  if (!fEventStat){
    fEventStat=new TH1D("hEventStat","Event statistics",5,0,5);
    fEventStat->GetXaxis()->SetBinLabel(1,"Before Phys. Sel.");
    fEventStat->GetXaxis()->SetBinLabel(2,"After Phys. Sel.");
    fEventStat->GetXaxis()->SetBinLabel(3,"After Cand. Sel.");
  }
  
  PostData(2,fEventStat);
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
    
    //ESD case
    if (man->GetInputEventHandler()->IsA()==AliESDInputHandler::Class()){
      if (!AliDielectronVarManager::GetESDpid()){
        
        if (AliDielectronMC::Instance()->HasMC()) {
          AliDielectronVarManager::InitESDpid();
        } else {
          AliDielectronVarManager::InitESDpid(1);
        }
      }
    }
    //AOD case
    if (man->GetInputEventHandler()->IsA()==AliAODInputHandler::Class()){
      if (!AliDielectronVarManager::GetAODpidUtil()){
        if (AliDielectronMC::Instance()->HasMC()) {
          AliDielectronVarManager::InitAODpidUtil();
        } else {
          AliDielectronVarManager::InitAODpidUtil(1);
        }
      }
    }
  }
  
  // Was event selected ?
  AliInputEventHandler* inputHandler = (AliInputEventHandler*) (man->GetInputEventHandler());
  UInt_t isSelected = AliVEvent::kAny;
  if( fSelectPhysics && inputHandler && inputHandler->GetEventSelection() ) {
    isSelected = inputHandler->IsEventSelected();
    isSelected&=fTriggerMask;
  }
  
  //Before physics selection
  fEventStat->Fill(0.);
  if (isSelected==0) {
    PostData(2,fEventStat);
    return;
  }
  //after physics selection
  fEventStat->Fill(1.);
  
  //bz for AliKF
  Double_t bz = InputEvent()->GetMagneticField();
  AliKFParticle::SetField( bz );
  
  fDielectron->Process(InputEvent());
  
  Bool_t hasCand = kFALSE;
  if(fStoreLikeSign) hasCand = (fDielectron->HasCandidates() || fDielectron->HasCandidatesLikeSign());
  else hasCand = (fDielectron->HasCandidates());
  
  if(hasCand){
    AliAODHandler *aodH=(AliAODHandler*)((AliAnalysisManager::GetAnalysisManager())->GetOutputEventHandler());
    AliAODEvent *aod = aodH->GetAOD();
    
    //replace the references of the legs with the AOD references
    TObjArray *obj = 0x0;
    for(Int_t i=0; i < 10; i++ ){
      obj = (TObjArray*)((*(fDielectron->GetPairArraysPointer()))->UncheckedAt(i));
      if(!obj) continue;
      for(int j=0;j<obj->GetEntriesFast();j++){
        AliAODTrack *leg1 = 0x0;
        AliAODTrack *leg2 = 0x0;
        AliDielectronPair *pairObj = (AliDielectronPair*)obj->UncheckedAt(j);
        Int_t id1 = ((AliVTrack*)pairObj->GetFirstDaughter())->GetID();
        Int_t id2 = ((AliVTrack*)pairObj->GetSecondDaughter())->GetID();
        
        for(Int_t it=0;it<aod->GetNumberOfTracks();it++){
          if(aod->GetTrack(it)->GetID() == id1) leg1 = aod->GetTrack(it);
          if(aod->GetTrack(it)->GetID() == id2) leg2 = aod->GetTrack(it);
        }
        if(!leg1 || !leg2) continue;
        
        if(man->GetInputEventHandler()->IsA()==AliAODInputHandler::Class()){
          leg1->ResetBit(kIsReferenced);
          leg1->SetUniqueID(0);
          leg2->ResetBit(kIsReferenced);
          leg2->SetUniqueID(0);
        }
        pairObj->SetRefFirstDaughter(leg1);
        pairObj->SetRefSecondDaughter(leg2);
      }
    }
    
    AliAODExtension *extDielectron = aodH->GetFilteredAOD("AliAOD.Dielectron.root");
    extDielectron->SelectEvent();
    //after candidate selection
    fEventStat->Fill(2.);
    
    //see if dielectron candidate branch exists, if not create is
    TTree *t=extDielectron->GetTree();

    if(!t->GetListOfBranches()->GetEntries() && man->GetInputEventHandler()->IsA()==AliAODInputHandler::Class())
      t->Branch(aod->GetList());
    
    if (!t->GetBranch("dielectrons")){
      t->Bronch("dielectrons","TObjArray",fDielectron->GetPairArraysPointer());
    }
    
    if(man->GetInputEventHandler()->IsA()==AliAODInputHandler::Class()) t->Fill();
  }
  
  PostData(1, const_cast<THashList*>(fDielectron->GetHistogramList()));
  PostData(2,fEventStat);
}

