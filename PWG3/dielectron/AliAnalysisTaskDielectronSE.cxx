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
//                      for Dielectron Analysis                          //
//                                                                       //
///////////////////////////////////////////////////////////////////////////

#include <TChain.h>
#include <TH1D.h>

#include <AliCFContainer.h>
#include <AliVEvent.h>
#include <AliInputEventHandler.h>
#include <AliESDInputHandler.h>
#include <AliAnalysisManager.h>
#include <AliAODInputHandler.h>

#include "AliDielectron.h"
#include "AliDielectronHistos.h"
#include "AliDielectronCF.h"
#include "AliAnalysisTaskDielectronSE.h"

ClassImp(AliAnalysisTaskDielectronSE)

//_________________________________________________________________________________
AliAnalysisTaskDielectronSE::AliAnalysisTaskDielectronSE() :
  AliAnalysisTaskSE(),
  fDielectron(0),
  fSelectPhysics(kFALSE),
  fTriggerMask(AliVEvent::kMB),
  fEventStat(0x0)
{
  //
  // Constructor
  //
}

//_________________________________________________________________________________
AliAnalysisTaskDielectronSE::AliAnalysisTaskDielectronSE(const char *name) :
  AliAnalysisTaskSE(name),
  fDielectron(0),
  fSelectPhysics(kFALSE),
  fTriggerMask(AliVEvent::kMB),
  fEventStat(0x0)
{
  //
  // Constructor
  //
  DefineInput(0,TChain::Class());
  DefineOutput(1, THashList::Class());
  DefineOutput(2, AliCFContainer::Class());
  DefineOutput(3, TH1D::Class());
}

//_________________________________________________________________________________
void AliAnalysisTaskDielectronSE::UserCreateOutputObjects()
{
  //
  // Initialise the framework objects
  //
  if (!fDielectron){
    AliError("No Dielectron framework object set !!!");
    return;
  }
  fDielectron->Init();
  if (fDielectron->GetHistogramList()){
    PostData(1, const_cast<THashList*>(fDielectron->GetHistogramList()));
  }
  if (fDielectron->GetCFManagerPair()){
    PostData(2, const_cast<AliCFContainer*>(fDielectron->GetCFManagerPair()->GetContainer()));
  }
  
  if (!fEventStat){
    fEventStat=new TH1D("hEventStat","Event statistics",5,0,5);
    fEventStat->GetXaxis()->SetBinLabel(1,"Before Phys. Sel.");
    fEventStat->GetXaxis()->SetBinLabel(2,"After Phys. Sel.");
  }
  
  PostData(3,fEventStat);
  
}

//_________________________________________________________________________________
void AliAnalysisTaskDielectronSE::UserExec(Option_t *)
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
    PostData(3,fEventStat);
    return;
  }
  //after physics selection
  fEventStat->Fill(1.);
  
  //bz for AliKF
  Double_t bz = InputEvent()->GetMagneticField();
  AliKFParticle::SetField( bz );
  
  fDielectron->Process(InputEvent());

  if (fDielectron->GetHistogramList()){
    PostData(1, const_cast<THashList*>(fDielectron->GetHistogramList()));
  }
  if (fDielectron->GetCFManagerPair()){
    PostData(2, const_cast<AliCFContainer*>(fDielectron->GetCFManagerPair()->GetContainer()));
  }
  PostData(3,fEventStat);
}

