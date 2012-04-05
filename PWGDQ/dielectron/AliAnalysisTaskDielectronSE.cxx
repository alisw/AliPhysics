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
#include <AliTriggerAnalysis.h>

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
  fTriggerOnV0AND(kFALSE),
  fRejectPileup(kFALSE),
  fTriggerAnalysis(0x0),
  fEventFilter(0x0),
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
  fTriggerOnV0AND(kFALSE),
  fRejectPileup(kFALSE),
  fTriggerAnalysis(0x0),
  fEventFilter(0x0),
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
  Int_t nbins=kNbinsEvent+2;
  if (!fEventStat){
    fEventStat=new TH1D("hEventStat","Event statistics",nbins,0,nbins);
    fEventStat->GetXaxis()->SetBinLabel(1,"Before Phys. Sel.");
    fEventStat->GetXaxis()->SetBinLabel(2,"After Phys. Sel.");

    //default names
    fEventStat->GetXaxis()->SetBinLabel(3,"Bin3 not used");
    fEventStat->GetXaxis()->SetBinLabel(4,"Bin4 not used");
    fEventStat->GetXaxis()->SetBinLabel(5,"Bin5 not used");

    if(fTriggerOnV0AND) fEventStat->GetXaxis()->SetBinLabel(3,"V0and triggers");
    if (fEventFilter) fEventStat->GetXaxis()->SetBinLabel(4,"After Event Filter");
    if (fRejectPileup) fEventStat->GetXaxis()->SetBinLabel(5,"After Pileup rejection");

    fEventStat->GetXaxis()->SetBinLabel((kNbinsEvent+1),Form("#splitline{1 candidate}{%s}",fDielectron->GetName()));
    fEventStat->GetXaxis()->SetBinLabel((kNbinsEvent+2),Form("#splitline{With >1 candidate}{%s}",fDielectron->GetName()));

  }

  if (!fTriggerAnalysis) fTriggerAnalysis=new AliTriggerAnalysis;
  fTriggerAnalysis->EnableHistograms();
  fTriggerAnalysis->SetAnalyzeMC(AliDielectronMC::Instance()->HasMC());
  
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
  Bool_t isESD=man->GetInputEventHandler()->IsA()==AliESDInputHandler::Class();
  Bool_t isAOD=man->GetInputEventHandler()->IsA()==AliAODInputHandler::Class();

  AliInputEventHandler* inputHandler = (AliInputEventHandler*) (man->GetInputEventHandler());
  if (!inputHandler) return;

  if ( inputHandler->GetPIDResponse() ){
    AliDielectronVarManager::SetPIDResponse( inputHandler->GetPIDResponse() );
  } else {
    AliFatal("This task needs the PID response attached to the input event handler!");
  }

  // Was event selected ?
  UInt_t isSelected = AliVEvent::kAny;
  if( fSelectPhysics && inputHandler){
  if((isESD && inputHandler->GetEventSelection()) || isAOD){
      isSelected = inputHandler->IsEventSelected();
      isSelected&=fTriggerMask;
      }
   }


  //Before physics selection
  fEventStat->Fill(kAllEvents);
  if (isSelected==0) {
    PostData(3,fEventStat);
    return;
  }
  //after physics selection
  fEventStat->Fill(kSelectedEvents);

  //V0and
  if(fTriggerOnV0AND){
  if(isESD){if (!fTriggerAnalysis->IsOfflineTriggerFired(static_cast<AliESDEvent*>(InputEvent()), AliTriggerAnalysis::kV0AND))
            return;}
  if(isAOD){if(!((static_cast<AliAODEvent*>(InputEvent()))->GetVZEROData()->GetV0ADecision() == AliVVZERO::kV0BB &&
            (static_cast<AliAODEvent*>(InputEvent()))->GetVZEROData()->GetV0CDecision() == AliVVZERO::kV0BB) )
            return;}
   }


  fEventStat->Fill(kV0andEvents);

  //Fill Event histograms before the event filter
  Double_t values[AliDielectronVarManager::kNMaxValues]={0};
  Double_t valuesMC[AliDielectronVarManager::kNMaxValues]={0};
  AliDielectronVarManager::Fill(InputEvent(),values);
  Bool_t hasMC=AliDielectronMC::Instance()->HasMC();
  if (hasMC) {
    if (AliDielectronMC::Instance()->ConnectMCEvent())
      AliDielectronVarManager::Fill(AliDielectronMC::Instance()->GetMCEvent(),valuesMC);
  }
  
  AliDielectronHistos *h=fDielectron->GetHistoManager();
  if (h){
    if (h->GetHistogramList()->FindObject("Event_noCuts"))
      h->FillClass("Event_noCuts",AliDielectronVarManager::kNMaxValues,values);
    if (hasMC && h->GetHistogramList()->FindObject("MCEvent_noCuts"))
      h->FillClass("Event_noCuts",AliDielectronVarManager::kNMaxValues,valuesMC);
  }
  
  //event filter
  if (fEventFilter) {
    if (!fEventFilter->IsSelected(InputEvent())) return;
  }
  fEventStat->Fill(kFilteredEvents);
  
  //pileup
  if (fRejectPileup){
    if (InputEvent()->IsPileupFromSPD(3,0.8,3.,2.,5.)) return;
  }
  fEventStat->Fill(kPileupEvents);
  
  //bz for AliKF
  Double_t bz = InputEvent()->GetMagneticField();
  AliKFParticle::SetField( bz );

  // make an artificial shift in the electron nsigma. Configured in the Config file
  AliDielectronPID::SetCorrVal((Double_t)InputEvent()->GetRunNumber());

  //
  // Actual data processing
  //
  fDielectron->Process(InputEvent());

  //statistics for number of selected candidates
  Int_t ncandidates=fDielectron->GetPairArray(1)->GetEntriesFast();
  if (ncandidates==1) fEventStat->Fill((kNbinsEvent));
  else if (ncandidates>1) fEventStat->Fill((kNbinsEvent+1));

  //Publish the data
  if (fDielectron->GetHistogramList()){
    PostData(1, const_cast<THashList*>(fDielectron->GetHistogramList()));
  }
  if (fDielectron->GetCFManagerPair()){
    PostData(2, const_cast<AliCFContainer*>(fDielectron->GetCFManagerPair()->GetContainer()));
  }
  PostData(3,fEventStat);
}

