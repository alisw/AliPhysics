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

#include <AliCFContainer.h>
#include <AliInputEventHandler.h>
#include <AliESDInputHandler.h>
#include <AliAnalysisManager.h>
#include <AliVEvent.h>

#include "AliDielectron.h"
#include "AliDielectronHistos.h"
#include "AliDielectronCF.h"
#include "AliDielectronMC.h"
#include "AliAnalysisTaskMultiDielectron.h"

ClassImp(AliAnalysisTaskMultiDielectron)

//_________________________________________________________________________________
AliAnalysisTaskMultiDielectron::AliAnalysisTaskMultiDielectron() :
  AliAnalysisTaskSE(),
  fListDielectron(),
  fListHistos(),
  fListCF(),
  fSelectPhysics(kFALSE)
{
  //
  // Constructor
  //
}

//_________________________________________________________________________________
AliAnalysisTaskMultiDielectron::AliAnalysisTaskMultiDielectron(const char *name) :
  AliAnalysisTaskSE(name),
  fListDielectron(),
  fListHistos(),
  fListCF(),
  fSelectPhysics(kFALSE)
{
  //
  // Constructor
  //
  DefineInput(0,TChain::Class());
  DefineOutput(1, TList::Class());
  DefineOutput(2, TList::Class());
  fListHistos.SetName("Dielectron_Histos_Multi");
  fListCF.SetName("Dielectron_CF_Multi");
}


//_________________________________________________________________________________
void AliAnalysisTaskMultiDielectron::UserCreateOutputObjects()
{
  //
  // Add all histogram manager histogram lists to the output TList
  //

  if (!fListHistos.IsEmpty()||!fListCF.IsEmpty()) return; //already initialised

  TIter nextDie(&fListDielectron);
  AliDielectron *die=0;
  while ( (die=static_cast<AliDielectron*>(nextDie())) ){
    die->Init();
    if (die->GetHistogramList()) fListHistos.Add(const_cast<THashList*>(die->GetHistogramList()));
    if (die->GetCFManagerPair()) fListCF.Add(const_cast<AliCFContainer*>(die->GetCFManagerPair()->GetContainer()));
  }
}

//_________________________________________________________________________________
void AliAnalysisTaskMultiDielectron::UserExec(Option_t *)
{
  //
  // Main loop. Called for every event
  //

  if (fListHistos.IsEmpty()&&fListCF.IsEmpty()) return;

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

  //Process event in all AliDielectron instances
  TIter nextDie(&fListDielectron);
  AliDielectron *die=0;
  while ( (die=static_cast<AliDielectron*>(nextDie())) ){
    die->Process(InputEvent());
  }
  
  PostData(1, &fListHistos);
  PostData(2, &fListCF);
}

//_________________________________________________________________________________
void AliAnalysisTaskMultiDielectron::FinishTaskOutput()
{
  //
  // Write debug tree
  //
  TIter nextDie(&fListDielectron);
  AliDielectron *die=0;
  while ( (die=static_cast<AliDielectron*>(nextDie())) ){
    die->SaveDebugTree();
  }
}

